#* =====第1步：corHMM模型比较、随机映射与逐事件定位=====

library(ape)
library(corHMM)
library(knitr)

#* =====读取命令行参数=====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || length(args) %% 2 != 0) {
  stop("参数必须以 --名称 值 成对传入")
}
arg_keys <- sub("^--", "", args[seq(1, length(args), by = 2)])
arg_values <- args[seq(2, length(args), by = 2)]
names(arg_values) <- arg_keys

required_keys <- c(
  "input-tree", "input-traits", "output-model", "output-simmap",
  "output-trans", "temp-dir", "id-column", "trait-column", "models",
  "criterion", "root-prior", "rate-categories", "node-states", "n-starts",
  "model-cores", "n-simulations", "simmap-cores", "max-attempt",
  "resolve-polytomies", "branch-floor-ratio", "scale-branch-lengths",
  "seed", "overwrite"
)
missing_keys <- setdiff(required_keys, names(arg_values))
if (length(missing_keys) > 0) {
  stop(sprintf("缺少参数：%s", paste(missing_keys, collapse = ", ")))
}

input_tree_path <- arg_values[["input-tree"]]
input_traits_path <- arg_values[["input-traits"]]
output_model_dir <- arg_values[["output-model"]]
output_simmap_dir <- arg_values[["output-simmap"]]
output_trans_dir <- arg_values[["output-trans"]]
temp_dir <- arg_values[["temp-dir"]]
id_column <- arg_values[["id-column"]]
trait_column <- arg_values[["trait-column"]]
model_names <- trimws(strsplit(arg_values[["models"]], ",", fixed = TRUE)[[1]])
selection_criterion <- toupper(arg_values[["criterion"]])
root_prior <- tolower(arg_values[["root-prior"]])
rate_categories <- as.integer(arg_values[["rate-categories"]])
node_states_method <- tolower(arg_values[["node-states"]])
n_starts <- as.integer(arg_values[["n-starts"]])
model_cores <- as.integer(arg_values[["model-cores"]])
n_simulations <- as.integer(arg_values[["n-simulations"]])
simmap_cores <- as.integer(arg_values[["simmap-cores"]])
max_attempt <- as.integer(arg_values[["max-attempt"]])
resolve_polytomies <- tolower(arg_values[["resolve-polytomies"]]) == "true"
branch_floor_ratio <- as.numeric(arg_values[["branch-floor-ratio"]])
scale_branch_lengths <- tolower(arg_values[["scale-branch-lengths"]]) == "true"
random_seed <- as.integer(arg_values[["seed"]])
overwrite <- tolower(arg_values[["overwrite"]]) == "true"

if (!all(model_names %in% c("ER", "SYM", "ARD"))) {
  stop("models仅允许ER、SYM和ARD")
}
if (length(unique(model_names)) != length(model_names)) {
  stop("models中存在重复模型")
}
if (!(selection_criterion %in% c("AIC", "BIC"))) {
  stop("criterion仅允许AIC或BIC")
}
if (!(root_prior %in% c("maddfitz", "yang", "flat"))) {
  stop("root-prior仅允许maddfitz、yang或flat")
}
if (rate_categories != 1) {
  stop("本示例仅支持rate_categories=1的标准Mk模型")
}
if (!(node_states_method %in% c("marginal", "joint"))) {
  stop("node-states仅允许marginal或joint")
}
if (any(is.na(c(n_starts, model_cores, n_simulations, simmap_cores,
                max_attempt, random_seed))) ||
    any(c(n_starts, model_cores, n_simulations, simmap_cores,
          max_attempt) < 1)) {
  stop("起点数、核心数、映射数和最大尝试次数必须为正整数")
}
if (!is.finite(branch_floor_ratio) || branch_floor_ratio <= 0) {
  stop("branch-floor-ratio必须为正数")
}

dir.create(output_model_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_simmap_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_trans_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

expected_outputs <- c(
  file.path(output_model_dir, "1-model-comparison.tsv"),
  file.path(output_model_dir, "1-best-model.tsv"),
  file.path(output_model_dir, "1-best-model-rate-matrix.tsv"),
  file.path(output_simmap_dir, "2-simmap-transition-summary.tsv"),
  file.path(output_simmap_dir, "2-simmap-edge-support.tsv"),
  file.path(output_simmap_dir, "2-simmap-maps.rds"),
  file.path(output_trans_dir, "3-transition-events.tsv"),
  file.path(output_trans_dir, "3-tree-edges.tsv"),
  file.path(output_trans_dir, "3-analysis.tree")
)
if (!overwrite && any(file.exists(expected_outputs))) {
  stop("输出已存在且overwrite=false")
}

set.seed(random_seed)

#* =====读取并校验输入=====
tree_input <- read.tree(input_tree_path)
dat1 <- read.delim(input_traits_path, stringsAsFactors = FALSE,
                   check.names = FALSE)
print(knitr::kable(head(dat1)))

if (!all(c(id_column, trait_column) %in% colnames(dat1))) {
  stop("性状表缺少配置指定的ID列或性状列")
}
if (ncol(dat1) != 2) {
  stop("示例性状表必须且只能包含ID和感兴趣性状两列")
}
dat1 <- dat1[, c(id_column, trait_column), drop = FALSE]
colnames(dat1) <- c("species", "state")
if (any(is.na(dat1$species) | dat1$species == "" |
        is.na(dat1$state) | dat1$state == "")) {
  stop("性状表的ID或性状存在缺失")
}
if (anyDuplicated(dat1$species) > 0) {
  stop("性状表存在重复ID")
}
if (anyDuplicated(tree_input$tip.label) > 0) {
  stop("树中存在重复末端ID")
}
if (!setequal(tree_input$tip.label, dat1$species)) {
  missing_traits <- setdiff(tree_input$tip.label, dat1$species)
  missing_tips <- setdiff(dat1$species, tree_input$tip.label)
  stop(sprintf(
    "树与性状ID不完全一致：缺性状%d个，缺树末端%d个",
    length(missing_traits), length(missing_tips)
  ))
}
if (!is.rooted(tree_input)) {
  stop("输入树必须有根；脚本不会擅自选择根")
}
if (is.null(tree_input$edge.length) || any(!is.finite(tree_input$edge.length)) ||
    any(tree_input$edge.length < 0)) {
  stop("输入树必须具有有限且非负的枝长")
}

dat1 <- dat1[match(tree_input$tip.label, dat1$species), , drop = FALSE]
state_levels <- sort(unique(dat1$state))
if (length(state_levels) < 2) {
  stop("至少需要两个性状状态")
}

cat(sprintf("输入：%d个末端，%d个状态\n", Ntip(tree_input), length(state_levels)))
print(table(dat1$state))

#* =====准备分析树=====
tree_analysis <- tree_input
if (!is.binary(tree_analysis)) {
  if (!resolve_polytomies) {
    stop("输入树含多分叉；请在配置中启用resolve_polytomies或先处理树")
  }
  tree_analysis <- multi2di(tree_analysis, random = FALSE)
}

edge_was_floored <- tree_analysis$edge.length <= 0
positive_lengths <- tree_analysis$edge.length[tree_analysis$edge.length > 0]
if (length(positive_lengths) == 0) {
  stop("树中没有正枝长")
}
branch_floor <- min(positive_lengths) * branch_floor_ratio
tree_analysis$edge.length[edge_was_floored] <- branch_floor
edge_length_input_units <- tree_analysis$edge.length
branch_scale_factor <- 1
if (scale_branch_lengths) {
  branch_scale_factor <- mean(tree_analysis$edge.length)
  tree_analysis$edge.length <- tree_analysis$edge.length / branch_scale_factor
}
if (!is.binary(tree_analysis) || any(tree_analysis$edge.length <= 0)) {
  stop("分析树必须完全二分且所有枝长为正")
}

write.tree(tree_analysis, file.path(output_trans_dir, "3-analysis.tree"))

n_tip <- Ntip(tree_analysis)
n_node_total <- n_tip + tree_analysis$Nnode
descendant_tips <- vector("list", n_node_total)
for (i in seq_len(n_tip)) descendant_tips[[i]] <- i
for (i in (n_tip + 1):n_node_total) descendant_tips[[i]] <- integer(0)
tree_postorder <- reorder(tree_analysis, "postorder")
for (e in seq_len(nrow(tree_postorder$edge))) {
  parent_node <- tree_postorder$edge[e, 1]
  child_node <- tree_postorder$edge[e, 2]
  descendant_tips[[parent_node]] <- c(
    descendant_tips[[parent_node]], descendant_tips[[child_node]]
  )
}

node_display <- ifelse(
  seq_len(n_node_total) <= n_tip,
  c(tree_analysis$tip.label, rep(NA_character_, tree_analysis$Nnode)),
  paste0("Node", seq_len(n_node_total))
)
edge_parent <- tree_analysis$edge[, 1]
edge_child <- tree_analysis$edge[, 2]
child_representative <- vapply(edge_child, function(nd) {
  tip_ids <- tree_analysis$tip.label[descendant_tips[[nd]]]
  paste(head(sort(tip_ids), 3), collapse = ";")
}, character(1))

dat2 <- data.frame(
  edge_id = seq_len(nrow(tree_analysis$edge)),
  parent_node = edge_parent,
  child_node = edge_child,
  parent_label = node_display[edge_parent],
  child_label = node_display[edge_child],
  child_is_tip = edge_child <= n_tip,
  n_descendant_tip = vapply(edge_child, function(nd) {
    length(descendant_tips[[nd]])
  }, integer(1)),
  descendant_tip_representative = child_representative,
  edge_length_analysis = tree_analysis$edge.length,
  edge_length_input_units = edge_length_input_units,
  edge_was_floored = edge_was_floored,
  branch_scale_factor = branch_scale_factor,
  stringsAsFactors = FALSE
)
write.table(dat2, file.path(output_trans_dir, "3-tree-edges.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

#* =====比较ER、SYM和ARD模型=====
corhmm_data <- data.frame(
  species = dat1$species,
  state = dat1$state,
  stringsAsFactors = FALSE
)
fit_objects <- list()
fit_records <- list()
n_state <- length(state_levels)

for (model_name in model_names) {
  cat(sprintf("拟合%s模型\n", model_name))
  fit <- try(
    corHMM(
      phy = tree_analysis,
      data = corhmm_data,
      rate.cat = rate_categories,
      model = model_name,
      root.p = root_prior,
      node.states = node_states_method,
      nstarts = n_starts,
      n.cores = model_cores
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    stop(sprintf("%s模型拟合失败：%s", model_name, as.character(fit)))
  }
  n_parameter <- switch(
    model_name,
    ER = 1L,
    SYM = as.integer(n_state * (n_state - 1) / 2),
    ARD = as.integer(n_state * (n_state - 1))
  )
  log_likelihood <- as.numeric(fit$loglik)
  aic <- as.numeric(fit$AIC)
  bic <- -2 * log_likelihood + n_parameter * log(n_tip)
  fit_objects[[model_name]] <- fit
  fit_records[[model_name]] <- data.frame(
    model = model_name,
    root_prior = root_prior,
    n_tip = n_tip,
    n_state = n_state,
    n_parameter = n_parameter,
    log_likelihood = log_likelihood,
    AIC = aic,
    BIC = bic,
    stringsAsFactors = FALSE
  )
}

dat3 <- do.call(rbind, fit_records)
rownames(dat3) <- NULL
criterion_values <- dat3[[selection_criterion]]
dat3$delta_AIC <- dat3$AIC - min(dat3$AIC)
dat3$delta_BIC <- dat3$BIC - min(dat3$BIC)
dat3$selection_criterion <- selection_criterion
dat3$selection_score <- criterion_values
dat3$rank <- rank(criterion_values, ties.method = "min")

best_candidates <- which(abs(criterion_values - min(criterion_values)) <= 1e-8)
best_row <- best_candidates[which.min(dat3$n_parameter[best_candidates])]
best_model_name <- dat3$model[best_row]
dat3$selected <- seq_len(nrow(dat3)) == best_row
dat3 <- dat3[order(dat3$rank, dat3$n_parameter), , drop = FALSE]

write.table(dat3, file.path(output_model_dir, "1-model-comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
dat4 <- dat3[dat3$selected, , drop = FALSE]
write.table(dat4, file.path(output_model_dir, "1-best-model.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

best_fit <- fit_objects[[best_model_name]]
q_matrix <- best_fit$solution
q_matrix[is.na(q_matrix)] <- 0
diag(q_matrix) <- 0
diag(q_matrix) <- -rowSums(q_matrix)
rownames(q_matrix) <- state_levels
colnames(q_matrix) <- state_levels
dat5 <- as.data.frame(as.table(q_matrix), stringsAsFactors = FALSE)
colnames(dat5) <- c("from_state", "to_state", "rate")
dat5$model <- best_model_name
dat5$root_prior <- root_prior
dat5 <- dat5[, c("model", "root_prior", "from_state", "to_state", "rate")]
write.table(dat5,
            file.path(output_model_dir, "1-best-model-rate-matrix.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("最优模型：%s（%s=%.6f）\n", best_model_name,
            selection_criterion, criterion_values[best_row]))

#* =====执行随机映射=====
# corHMM 2.8在makeSimmap(root.p="flat")内部错误地把矩阵维度向量传给rep。
# 传入标量1会对所有根状态施加相同权重，与flat先验在归一化后完全等价。
simmap_root_prior <- if (root_prior == "flat") 1 else root_prior
simmap_result <- makeSimmap(
  tree = tree_analysis,
  data = corhmm_data,
  model = q_matrix,
  rate.cat = rate_categories,
  root.p = simmap_root_prior,
  nSim = n_simulations,
  nCores = simmap_cores,
  max.attempt = max_attempt
)
if (length(simmap_result) != n_simulations) {
  stop("随机映射返回数量与配置不一致")
}
saveRDS(simmap_result, file.path(output_simmap_dir, "2-simmap-maps.rds"))

#* =====提取每次模拟的逐事件位置=====
event_records <- list()
event_record_index <- 0L

for (simulation_id in seq_len(n_simulations)) {
  one_map <- simmap_result[[simulation_id]]
  if (length(one_map$maps) != nrow(tree_analysis$edge)) {
    stop(sprintf("第%d次映射的枝数量不一致", simulation_id))
  }

  # 校验内部节点两侧状态连续，防止误读maps中的状态顺序。
  for (node_id in (n_tip + 1):n_node_total) {
    incoming_edge <- which(edge_child == node_id)
    outgoing_edges <- which(edge_parent == node_id)
    if (length(incoming_edge) == 1 && length(outgoing_edges) > 0) {
      incoming_state <- tail(names(one_map$maps[[incoming_edge]]), 1)
      outgoing_states <- vapply(outgoing_edges, function(e) {
        head(names(one_map$maps[[e]]), 1)
      }, character(1))
      if (any(outgoing_states != incoming_state)) {
        stop(sprintf("第%d次映射在节点%d处状态不连续", simulation_id, node_id))
      }
    }
  }

  for (edge_id in seq_len(nrow(tree_analysis$edge))) {
    segments <- one_map$maps[[edge_id]]
    if (abs(sum(segments) - tree_analysis$edge.length[edge_id]) > 1e-8) {
      stop(sprintf("第%d次映射第%d条枝的分段长度不守恒",
                   simulation_id, edge_id))
    }
    if (length(segments) < 2) next

    state_path <- names(segments)
    event_positions <- cumsum(as.numeric(segments))[-length(segments)]
    for (event_order in seq_along(event_positions)) {
      event_record_index <- event_record_index + 1L
      position_analysis <- event_positions[event_order]
      relative_position <- position_analysis / tree_analysis$edge.length[edge_id]
      event_records[[event_record_index]] <- data.frame(
        simulation_id = simulation_id,
        event_id = sprintf("sim%04d_edge%04d_event%02d",
                           simulation_id, edge_id, event_order),
        edge_id = edge_id,
        parent_node = edge_parent[edge_id],
        child_node = edge_child[edge_id],
        parent_label = node_display[edge_parent[edge_id]],
        child_label = node_display[edge_child[edge_id]],
        event_order_on_edge = event_order,
        from_state = state_path[event_order],
        to_state = state_path[event_order + 1L],
        edge_length_analysis = tree_analysis$edge.length[edge_id],
        distance_from_parent_analysis = position_analysis,
        distance_to_child_analysis = tree_analysis$edge.length[edge_id] -
          position_analysis,
        relative_position_from_parent = relative_position,
        edge_length_input_units = edge_length_input_units[edge_id],
        distance_from_parent_input_units = position_analysis * branch_scale_factor,
        distance_to_child_input_units =
          (tree_analysis$edge.length[edge_id] - position_analysis) *
          branch_scale_factor,
        n_descendant_tip = dat2$n_descendant_tip[edge_id],
        descendant_tip_representative =
          dat2$descendant_tip_representative[edge_id],
        edge_was_floored = edge_was_floored[edge_id],
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(event_records) == 0) {
  stop("所有随机映射均未产生状态转移，无法输出逐事件结果")
}
dat6 <- do.call(rbind, event_records)
rownames(dat6) <- NULL
write.table(dat6, file.path(output_trans_dir, "3-transition-events.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

#* =====汇总方向性转移次数=====
dat6$transition <- paste(dat6$from_state, dat6$to_state, sep = "->")
transition_types <- sort(unique(dat6$transition))
transition_summary_records <- list()

for (transition_type in transition_types) {
  event_subset <- dat6[dat6$transition == transition_type, , drop = FALSE]
  event_counts <- tabulate(event_subset$simulation_id, nbins = n_simulations)
  transition_parts <- strsplit(transition_type, "->", fixed = TRUE)[[1]]
  transition_summary_records[[transition_type]] <- data.frame(
    model = best_model_name,
    root_prior = root_prior,
    from_state = transition_parts[1],
    to_state = transition_parts[2],
    mean_n_event = mean(event_counts),
    median_n_event = median(event_counts),
    sd_n_event = sd(event_counts),
    lower_n_event = as.numeric(quantile(event_counts, 0.025)),
    upper_n_event = as.numeric(quantile(event_counts, 0.975)),
    proportion_nonzero = mean(event_counts > 0),
    n_simulations = n_simulations,
    stringsAsFactors = FALSE
  )
}

dat7 <- do.call(rbind, transition_summary_records)
rownames(dat7) <- NULL
write.table(dat7,
            file.path(output_simmap_dir, "2-simmap-transition-summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

#* =====汇总逐枝方向性支持率=====
edge_transition_keys <- unique(dat6[, c("edge_id", "from_state", "to_state")])
edge_support_records <- list()

for (i in seq_len(nrow(edge_transition_keys))) {
  edge_id <- edge_transition_keys$edge_id[i]
  from_state <- edge_transition_keys$from_state[i]
  to_state <- edge_transition_keys$to_state[i]
  event_subset <- dat6[
    dat6$edge_id == edge_id & dat6$from_state == from_state &
      dat6$to_state == to_state,
    , drop = FALSE
  ]
  event_counts <- tabulate(event_subset$simulation_id, nbins = n_simulations)
  edge_support_records[[i]] <- data.frame(
    edge_id = edge_id,
    parent_node = dat2$parent_node[edge_id],
    child_node = dat2$child_node[edge_id],
    parent_label = dat2$parent_label[edge_id],
    child_label = dat2$child_label[edge_id],
    from_state = from_state,
    to_state = to_state,
    n_simulations_with_event = sum(event_counts > 0),
    edge_transition_support = mean(event_counts > 0),
    mean_n_event = mean(event_counts),
    n_descendant_tip = dat2$n_descendant_tip[edge_id],
    descendant_tip_representative =
      dat2$descendant_tip_representative[edge_id],
    edge_was_floored = dat2$edge_was_floored[edge_id],
    stringsAsFactors = FALSE
  )
}

dat8 <- do.call(rbind, edge_support_records)
dat8 <- dat8[order(-dat8$edge_transition_support, dat8$edge_id,
                   dat8$from_state, dat8$to_state), , drop = FALSE]
rownames(dat8) <- NULL
write.table(dat8, file.path(output_simmap_dir, "2-simmap-edge-support.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf(
  "完成：%d个模型，%d次随机映射，%d个逐事件，%d条逐枝方向记录\n",
  nrow(dat3), n_simulations, nrow(dat6), nrow(dat8)
))
print(dat4)
