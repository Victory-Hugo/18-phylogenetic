# --- 基础依赖 ---
library(jsonlite)
library(ape)
library(treeio)
library(ggtree)
library(dplyr)
library(ggplot2)
library(colorspace)
library(ggtreeExtra)

`%||%` <- function(x, y) if (!is.null(x)) x else y

make_safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) {
    x <- "item"
  }
  x
}

make_unique_tip_label <- function(label, existing) {
  if (!label %in% existing) {
    return(label)
  }
  suffix <- 1L
  candidate <- paste0(label, "_placement", suffix)
  while (candidate %in% existing) {
    suffix <- suffix + 1L
    candidate <- paste0(label, "_placement", suffix)
  }
  candidate
}

construct_augmented_tree <- function(base_phylo, placements_tbl) {
  if (is.null(placements_tbl) || nrow(placements_tbl) == 0) {
    return(list(phylo = base_phylo, mapping = data.frame(original = character(), final = character(), stringsAsFactors = FALSE)))
  }

  placements_tbl <- placements_tbl %>%
    mutate(
      distal_length = ifelse(is.na(distal_length), 0, distal_length),
      pendant_length = ifelse(is.na(pendant_length), 0, pendant_length)
    )

  edge_num_attr <- attr(base_phylo, "edgeNum")
  edge_map_idx <- match(edge_num_attr$node, base_phylo$edge[, 2])
  edge_num_for_edge <- rep(NA_real_, nrow(base_phylo$edge))
  valid_idx <- !is.na(edge_map_idx)
  edge_num_for_edge[edge_map_idx[valid_idx]] <- edge_num_attr$edgeNum[valid_idx]

  placements_tbl <- placements_tbl %>% filter(edge_num %in% edge_num_for_edge)
  if (nrow(placements_tbl) == 0) {
    return(list(phylo = base_phylo, mapping = data.frame(original = character(), final = character(), stringsAsFactors = FALSE)))
  }

  placements_split <- split(placements_tbl, placements_tbl$edge_num)

  new_tip_labels <- base_phylo$tip.label
  tip_node_ids <- seq_len(length(base_phylo$tip.label))
  new_edges <- list()
  new_lengths <- numeric()
  has_lengths <- !is.null(base_phylo$edge.length)
  next_node_index <- max(base_phylo$edge) + 1
  new_internal_count <- 0L
  mapping_entries <- list()

  add_edge <- function(parent, child, length_value) {
    new_edges[[length(new_edges) + 1]] <<- c(parent, child)
    if (has_lengths) {
      new_lengths <<- c(new_lengths, length_value)
    }
  }

  for (i in seq_len(nrow(base_phylo$edge))) {
    parent <- base_phylo$edge[i, 1]
    child <- base_phylo$edge[i, 2]
    base_length <- if (has_lengths) base_phylo$edge.length[i] else NA_real_
    edge_num_value <- edge_num_for_edge[i]
    placement_list <- placements_split[[as.character(edge_num_value)]]

    if (!is.null(placement_list) && nrow(placement_list) > 0) {
      placement_list <- placement_list %>% arrange(distal_length)
      current_child <- child
      current_distance <- 0
      for (j in seq_len(nrow(placement_list))) {
        p_row <- placement_list[j, ]
        d <- p_row$distal_length
        if (has_lengths && !is.na(base_length)) {
          d <- max(min(d, base_length), current_distance)
        } else {
          d <- max(d, current_distance)
        }
        segment_to_child <- d - current_distance
        if (segment_to_child < 0) segment_to_child <- 0
        new_internal <- next_node_index
        next_node_index <- next_node_index + 1L
        new_internal_count <- new_internal_count + 1L
        final_label <- make_unique_tip_label(p_row$name, new_tip_labels)
        new_tip_labels <- c(new_tip_labels, final_label)
        new_tip_node <- next_node_index
        next_node_index <- next_node_index + 1L
        tip_node_ids <- c(tip_node_ids, new_tip_node)

        add_edge(new_internal, current_child, segment_to_child)
        add_edge(new_internal, new_tip_node, p_row$pendant_length)

        mapping_entries[[length(mapping_entries) + 1]] <- c(original = p_row$name, final = final_label)
        current_child <- new_internal
        current_distance <- d
      }
      parent_segment <- if (has_lengths && !is.na(base_length)) max(base_length - current_distance, 0) else NA_real_
      add_edge(parent, current_child, parent_segment)
    } else {
      add_edge(parent, child, base_length)
    }
  }

  new_edge_matrix <- do.call(rbind, new_edges)
  unique_nodes <- sort(unique(c(new_edge_matrix)))
  tip_mapping <- setNames(seq_along(tip_node_ids), as.character(tip_node_ids))
  internal_nodes_ids <- setdiff(unique_nodes, tip_node_ids)
  internal_mapping <- setNames(length(tip_node_ids) + seq_along(internal_nodes_ids), as.character(internal_nodes_ids))
  node_mapping <- c(tip_mapping, internal_mapping)
  remapped_values <- node_mapping[match(new_edge_matrix, names(node_mapping))]
  if (any(is.na(remapped_values))) {
    stop("Failed to remap node indices when constructing augmented tree.")
  }
  edge_renumbered <- matrix(as.integer(remapped_values), ncol = 2)

  result_phylo <- list(
    edge = edge_renumbered,
    edge.length = if (has_lengths) new_lengths else NULL,
    tip.label = new_tip_labels,
    Nnode = length(internal_nodes_ids),
    root.edge = base_phylo$root.edge
  )
  class(result_phylo) <- "phylo"
  result_phylo <- ape::reorder.phylo(result_phylo)

  mapping_df <- if (length(mapping_entries) > 0) {
    as.data.frame(do.call(rbind, mapping_entries), stringsAsFactors = FALSE)
  } else {
    data.frame(original = character(), final = character(), stringsAsFactors = FALSE)
  }

  list(phylo = result_phylo, mapping = mapping_df)
}

# 从命令行参数读取输入
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("用法: Rscript 2-可视化.R <jplace文件路径> [分类文件路径] [输出目录]\n")
  cat("示例: Rscript 2-可视化.R data.jplace classification.csv output/\n")
  quit(status = 1)
}

jplace_file <- args[1]
if (!file.exists(jplace_file)) {
  cat("错误：jplace文件不存在:", jplace_file, "\n")
  quit(status = 1)
}

# 分类文件（可选）
classification_file <- if (length(args) >= 2 && nzchar(args[2])) args[2] else ""

# 输出目录（可选）
if (length(args) >= 3 && nzchar(args[3])) {
  output_dir <- args[3]
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
} else {
  output_dir <- dirname(jplace_file)
}

# 序列目标（默认为空）
sequence_targets <- character(0)

cat("输入文件:", jplace_file, "\n")
cat("输出目录:", output_dir, "\n")
if (nzchar(classification_file)) {
  cat("分类文件:", classification_file, "\n")
}

jplace_raw    <- jsonlite::read_json(jplace_file, simplifyVector = FALSE)
jplace_simple <- jsonlite::fromJSON(jplace_file)

to_scalar <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.atomic(x) && length(x) == 1) return(as.character(x))
  jsonlite::toJSON(x, auto_unbox = TRUE)
}

# 拆解 jplace 到多个 CSV 文件
overview_entries <- list(
  version = jplace_raw$version %||% NA_character_,
  tree    = jplace_raw$tree %||% NA_character_
)
overview_df <- data.frame(
  key = names(overview_entries),
  value = vapply(overview_entries, to_scalar, character(1)),
  stringsAsFactors = FALSE
)
write.csv(overview_df, file = file.path(output_dir, "jplace_overview.csv"), row.names = FALSE)

metadata_list <- jplace_raw$metadata
if (!is.null(metadata_list) && length(metadata_list) > 0) {
  metadata_df <- data.frame(
    key = names(metadata_list),
    value = vapply(metadata_list, to_scalar, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(metadata_df, file = file.path(output_dir, "jplace_metadata.csv"), row.names = FALSE)
}

fields_vec <- jplace_raw$fields %||% character()
if (length(fields_vec) > 0) {
  fields_df <- data.frame(
    field_index = seq_along(fields_vec),
    field_name = fields_vec,
    stringsAsFactors = FALSE
  )
  write.csv(fields_df, file = file.path(output_dir, "jplace_fields.csv"), row.names = FALSE)
}

placements_list <- jplace_raw$placements %||% list()

placement_rows <- lapply(seq_along(placements_list), function(pid) {
  entry <- placements_list[[pid]]
  p <- entry$p
  if (is.null(p) || length(p) == 0) return(NULL)
  rows <- lapply(p, function(row) {
    if (is.null(row)) rep(NA_real_, length(fields_vec)) else as.numeric(row)
  })
  candidate_df <- as.data.frame(do.call(rbind, rows))
  colnames(candidate_df) <- if (length(fields_vec) >= ncol(candidate_df)) {
    fields_vec[seq_len(ncol(candidate_df))]
  } else {
    paste0("field_", seq_len(ncol(candidate_df)))
  }
  candidate_df$placement_id   <- pid
  candidate_df$candidate_rank <- seq_len(nrow(candidate_df))
  candidate_df
})

placement_rows <- placement_rows[lengths(placement_rows) > 0]
placement_candidates <- if (length(placement_rows) > 0) {
  dplyr::bind_rows(placement_rows) %>% dplyr::relocate(placement_id, candidate_rank)
} else {
  tibble::tibble()
}

if (nrow(placement_candidates) > 0) {
  write.csv(placement_candidates, file = file.path(output_dir, "jplace_placement_candidates.csv"), row.names = FALSE)
}

name_rows <- lapply(seq_along(placements_list), function(pid) {
  entry <- placements_list[[pid]]
  nm <- entry$nm
  n  <- entry$n
  if (!is.null(nm) && length(nm) > 0) {
    nm_df <- lapply(seq_along(nm), function(idx) {
      pair <- nm[[idx]]
      data.frame(
        placement_id = pid, name_rank = idx,
        name = to_scalar(pair[[1]]),
        multiplicity = if (length(pair) >= 2) as.numeric(pair[[2]]) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(nm_df)
  } else if (!is.null(n) && length(n) > 0) {
    n_df <- lapply(seq_along(n), function(idx) {
      data.frame(placement_id = pid, name_rank = idx, name = to_scalar(n[[idx]]), 
                 multiplicity = NA_real_, stringsAsFactors = FALSE)
    })
    dplyr::bind_rows(n_df)
  }
})

name_rows <- name_rows[lengths(name_rows) > 0]
if (length(name_rows) > 0) {
  placement_names <- dplyr::bind_rows(name_rows)
  write.csv(placement_names, file = file.path(output_dir, "jplace_placement_names.csv"), row.names = FALSE)
}

extra_rows <- lapply(seq_along(placements_list), function(pid) {
  entry <- placements_list[[pid]]
  extra_keys <- setdiff(names(entry), c("p", "n", "nm"))
  if (length(extra_keys) == 0) return(NULL)
  data.frame(
    placement_id = pid, key = extra_keys,
    value = vapply(extra_keys, function(k) to_scalar(entry[[k]]), character(1)),
    stringsAsFactors = FALSE
  )
})

extra_rows <- extra_rows[lengths(extra_rows) > 0]
if (length(extra_rows) > 0) {
  placement_extras <- dplyr::bind_rows(extra_rows)
  write.csv(placement_extras, file = file.path(output_dir, "jplace_placement_extras.csv"), row.names = FALSE)
}

# 读取 jplace 并构造 phylo
# 优先使用 treeio::read.jplace；若失败则从原始 JSON 降级处理
jobj <- tryCatch(treeio::read.jplace(jplace_file), error = function(e) NULL)
if (!is.null(jobj)) {
  phylo <- jobj@phylo
  placements <- tryCatch({
    suppressWarnings(as_tibble(jobj, what = "placement"))
  }, error = function(e) NULL)
  
  if (is.null(placements) || !"like_weight_ratio" %in% names(placements)) {
    placements <- getFromNamespace("extract.placement", "treeio")(jplace_simple, phylo)
    placements <- as_tibble(placements)
  }
} else {
  newick_clean <- gsub("\\{[^}]*\\}", "", jplace_simple$tree)
  phylo <- ape::read.tree(text = newick_clean)
  placements <- getFromNamespace("extract.placement", "treeio")(jplace_simple, phylo)
  placements <- as_tibble(placements)
}
write.csv(
  placements,
  file = file.path(output_dir, "df_pplacer_likelihood.csv"),
  row.names = FALSE
)

# 过滤：每条序列保留 like_weight_ratio 最大的放置
filtered_placements <- placements %>%
  group_by(name) %>%
  filter(like_weight_ratio == max(like_weight_ratio, na.rm = TRUE)) %>%
  ungroup()

# 直方图：like_weight_ratio 分布
dat <- bind_rows(
  transmute(placements, group = "Before",  likelihood_weight_ratio = like_weight_ratio),
  transmute(filtered_placements, group = "Filtered", likelihood_weight_ratio = like_weight_ratio)
)

p1 <- ggplot(dat, aes(x = likelihood_weight_ratio, fill = group)) +
  geom_histogram(position = "dodge", color = "white", linewidth = 0.5, binwidth = 0.1) +
  scale_fill_manual(values = c("Before" = "#f08c8c", "Filtered" = "#8dc5fe")) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  labs(x = "Likelihood Weight Ratio", y = "Number of Placements") +
  theme_bw()
ggplot2::ggsave(filename = file.path(output_dir, "likelihood_weight_ratio_hist.png"), plot = p1, width = 8, height = 5, dpi = 300)

# treedata：把节点放置数合并入树
tree_td <- as.treedata(phylo)
td_tbl  <- as_tibble(tree_td)
n_tip   <- length(phylo$tip.label)

node_counts <- filtered_placements %>%
  count(node, name = "nplace")

placement_labels <- filtered_placements %>%
  group_by(node) %>%
  summarise(placement_names = paste(sort(unique(name)), collapse = ", "), .groups = "drop")

post_prob_summary <- placements %>%
  group_by(node) %>%
  summarise(post_prob_max = if (all(is.na(post_prob))) NA_real_ else max(post_prob, na.rm = TRUE), .groups = "drop") %>%
  mutate(post_prob_max = ifelse(is.infinite(post_prob_max), NA_real_, post_prob_max))

# 合并数据到 treedata
tree_td_data <- td_tbl %>%
  left_join(node_counts, by = "node") %>%
  left_join(placement_labels, by = "node") %>%
  left_join(post_prob_summary, by = "node") %>%
  mutate(
    nplace = ifelse(is.na(.data$nplace), 0L, .data$nplace),
    label_original = dplyr::coalesce(.data$label, ""),
    label_display = label_original,
    post_prob_label = ifelse(!is.na(post_prob_max), sprintf("%.3f", post_prob_max), NA_character_),
    label = label_display
  ) %>%
  tibble::as_tibble()

tree_td@data <- tree_td_data

# 分类信息处理（可选）
classification_for_tree <- NULL
if (nzchar(classification_file) && file.exists(classification_file)) {
  classification_tbl <- tryCatch(
    read.csv(classification_file, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      message("读取分类信息失败：", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(classification_tbl) && nrow(classification_tbl) > 0) {
    label_candidates <- intersect(c("label", "tip", "tip_label", "tip.label", "taxon", "name"), names(classification_tbl))
    group_candidates <- intersect(c("group", "Group", "taxon_group", "category", "clade"), names(classification_tbl))
    
    if (length(label_candidates) > 0 && length(group_candidates) > 0) {
      label_col <- label_candidates[[1]]
      group_col <- group_candidates[[1]]
      classification_for_tree <- tree_td_data %>%
        filter(.data$node <= n_tip) %>%
        select(label_display, label_original) %>%
        left_join(classification_tbl, by = setNames(label_col, "label_original")) %>%
        mutate(group = .data[[group_col]], group = ifelse(is.na(group) | group == "", NA_character_, as.character(group)))
      
      if (all(is.na(classification_for_tree$group))) {
        classification_for_tree <- NULL
        message("分类文件已读取，但未找到有效的 group 信息。")
      }
    } else {
      message("分类文件缺少 label/group 字段。")
    }
  }
}

# 圆形树可视化：节点颜色/大小 = 放置数
max_n <- max(tree_td@data$nplace %||% 0L, na.rm = TRUE)
max_n <- if (is.finite(max_n) && max_n > 0) max_n else 1
color_gradient <- c("#BFBFBF", colorspace::sequential_hcl(6, "Sunset"))
branch_modes <- list(ignore_branch_length = "none", with_branch_length = "branch.length")
generated_main_files <- character(0)
grouped_files <- character(0)

for (mode_name in names(branch_modes)) {
  branch_setting <- branch_modes[[mode_name]]
  p2 <- ggtree(tree_td, layout = "circular", branch.length = branch_setting, lwd = 1, aes(color = nplace, size = nplace)) +
    scale_color_gradientn(colors = color_gradient, limits = c(0, max_n), na.value = "#D9D9D9") +
    scale_size_continuous(range = c(0.2, 2.6), limits = c(0, max_n), guide = guide_legend(reverse = TRUE)) +
    geom_tiplab2(aes(label = label_display), offset = 0.3, size = 2.2, show.legend = FALSE) +
    geom_text2(aes(subset = !is.na(post_prob_label), label = post_prob_label), size = 1.8, show.legend = FALSE) +
    theme(legend.position = "right", legend.box = "horizontal", legend.box.just = "left")
  
  main_file <- file.path(output_dir, paste0("placement_circular_tree_", mode_name, ".png"))
  ggplot2::ggsave(filename = main_file, plot = p2, width = 8, height = 8, dpi = 300)
  generated_main_files <- c(generated_main_files, main_file)
  
  if (!is.null(classification_for_tree)) {
    classification_layer <- classification_for_tree %>% mutate(label = label_display)
    group_levels <- classification_layer$group %>% unique() %>% setdiff(NA_character_) %>% .[. != ""]
    
    if (length(group_levels) > 0) {
      palette_groups <- colorspace::qualitative_hcl(length(group_levels), palette = "Dark 3")
      p2_groups <- p2 +
        ggtreeExtra::geom_fruit(data = classification_layer, geom = geom_tile, 
                                 mapping = aes(fill = factor(group, levels = group_levels)), 
                                 width = 0.6, offset = 0.15, inherit.aes = FALSE) +
        scale_fill_manual(name = "Taxon Group", values = palette_groups, 
                          guide = guide_legend(keywidth = 0.6, keyheight = 0.6, order = 3), na.translate = FALSE)
      
      grouped_file <- file.path(output_dir, paste0("placement_circular_tree_groups_", mode_name, ".png"))
      ggplot2::ggsave(filename = grouped_file, plot = p2_groups, width = 9, height = 8.5, dpi = 300)
      grouped_files <- c(grouped_files, grouped_file)
    }
  }
}

# 清理旧版本文件
legacy_main <- file.path(output_dir, "placement_circular_tree.png")
if (file.exists(legacy_main) && !(legacy_main %in% generated_main_files)) {
  invisible(file.remove(legacy_main))
}
legacy_grouped <- file.path(output_dir, "placement_circular_tree_groups.png")
if (file.exists(legacy_grouped) && !(legacy_grouped %in% grouped_files)) {
  invisible(file.remove(legacy_grouped))
}

# 指定序列放置不确定性展示
sequence_targets_input <- if (length(sequence_targets) > 0) {
  sequence_targets[!is.na(sequence_targets) & nzchar(sequence_targets)]
} else {
  unique(filtered_placements$name)
}
sequence_targets_input <- intersect(sequence_targets_input, unique(placements$name))

existing_sequence_plots <- list.files(output_dir, pattern = "^placement_sequence_.*\\.png$", full.names = TRUE)
if (length(existing_sequence_plots) > 0) invisible(file.remove(existing_sequence_plots))

if (length(sequence_targets_input) > 0) {
  for (seq_name in sequence_targets_input) {
    seq_data <- placements %>% filter(name == seq_name)
    if (nrow(seq_data) == 0) next
    
    seq_summary <- seq_data %>%
      group_by(node) %>%
      summarise(like_weight_ratio = max(like_weight_ratio, na.rm = TRUE), 
                post_prob = max(post_prob, na.rm = TRUE), .groups = "drop") %>%
      mutate(across(where(is.numeric), ~ ifelse(is.infinite(.), NA_real_, .)))
    
    seq_scale_max <- max(seq_summary$like_weight_ratio, na.rm = TRUE)
    seq_scale_max <- if (is.finite(seq_scale_max) && seq_scale_max > 0) seq_scale_max else 1
    
    for (mode_name in names(branch_modes)) {
      branch_setting <- branch_modes[[mode_name]]
      seq_tree_data <- tree_td_data %>%
        left_join(seq_summary, by = "node") %>%
        mutate(post_prob_seq = ifelse(!is.na(post_prob), post_prob, post_prob_max),
               post_prob_seq = ifelse(is.infinite(post_prob_seq), NA_real_, post_prob_seq),
               post_prob_label_seq = ifelse(!is.na(post_prob_seq), sprintf("%.3f", post_prob_seq), NA_character_))
      
      seq_tree <- tree_td
      seq_tree@data <- seq_tree_data
      
      p_seq <- ggtree(seq_tree, layout = "circular", branch.length = branch_setting, aes(color = like_weight_ratio)) +
        scale_color_gradientn(colors = color_gradient, limits = c(0, seq_scale_max), na.value = "#D9D9D9") +
        geom_point2(aes(subset = !is.na(like_weight_ratio)), size = 1.4) +
        geom_tiplab2(aes(label = label_display), offset = 0.25, size = 2, show.legend = FALSE) +
        geom_text2(aes(subset = !is.na(post_prob_label_seq), label = post_prob_label_seq), size = 1.6, show.legend = FALSE) +
        labs(title = paste0("Placement weights for ", seq_name), color = "Like weight ratio") +
        theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"))
      
      seq_file <- file.path(output_dir, paste0("placement_sequence_", make_safe_filename(seq_name), "_", mode_name, ".png"))
      ggplot2::ggsave(filename = seq_file, plot = p_seq, width = 7, height = 7, dpi = 300)
    }
  }
}

# 节点分支权重可视化
branch_summary <- placements %>%
  group_by(node) %>%
  summarise(branch_lwr = max(like_weight_ratio, na.rm = TRUE), 
            branch_post_prob = max(post_prob, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ ifelse(is.infinite(.), NA_real_, .)))

branch_tree_data <- tree_td_data %>%
  left_join(branch_summary, by = "node") %>%
  mutate(branch_post_label = ifelse(!is.na(branch_post_prob), sprintf("%.3f", branch_post_prob), NA_character_))

branch_scale_max <- max(branch_summary$branch_lwr, na.rm = TRUE)
branch_scale_max <- if (is.finite(branch_scale_max) && branch_scale_max > 0) branch_scale_max else 1

generated_branch_files <- character(0)
for (mode_name in names(branch_modes)) {
  branch_setting <- branch_modes[[mode_name]]
  branch_tree <- tree_td
  branch_tree@data <- branch_tree_data
  
  p_branch <- ggtree(branch_tree, layout = "circular", branch.length = branch_setting, aes(color = branch_lwr, size = branch_lwr)) +
    scale_color_gradientn(colors = color_gradient, limits = c(0, branch_scale_max), na.value = "#D9D9D9") +
    scale_size_continuous(range = c(0.2, 2.6), limits = c(0, branch_scale_max), guide = guide_legend(reverse = TRUE)) +
    geom_tiplab2(aes(label = label_display), offset = 0.3, size = 2.2, show.legend = FALSE) +
    geom_text2(aes(subset = !is.na(branch_post_label), label = branch_post_label), size = 1.8, show.legend = FALSE) +
    labs(color = "Max LWR", size = "Max LWR") +
    theme(legend.position = "right", legend.box = "horizontal", legend.box.just = "left")
  
  branch_file <- file.path(output_dir, paste0("placement_branch_weights_", mode_name, ".png"))
  ggplot2::ggsave(filename = branch_file, plot = p_branch, width = 8, height = 8, dpi = 300)
  generated_branch_files <- c(generated_branch_files, branch_file)
}

legacy_branch <- file.path(output_dir, "placement_branch_weights.png")
if (file.exists(legacy_branch) && !(legacy_branch %in% generated_branch_files)) {
  invisible(file.remove(legacy_branch))
}

# 构造增强树（包含放置信息）
best_insertions <- filtered_placements %>%
  arrange(desc(like_weight_ratio), desc(post_prob), desc(marginal_like)) %>%
  group_by(name) %>%
  slice_head(n = 1) %>%
  ungroup()

augmented_result <- construct_augmented_tree(phylo, best_insertions)
augmented_phylo <- augmented_result$phylo
ape::write.tree(augmented_phylo, file = file.path(output_dir, "placement_augmented_tree.nwk"))

mapping_df <- augmented_result$mapping
mapping_file <- file.path(output_dir, "placement_augmented_mapping.csv")
if (nrow(mapping_df) > 0) {
  write.csv(mapping_df, file = mapping_file, row.names = FALSE)
} else if (file.exists(mapping_file)) {
  invisible(file.remove(mapping_file))
}
