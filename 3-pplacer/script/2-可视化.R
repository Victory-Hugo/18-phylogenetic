# --- 基础依赖（最小集合） ---
library(jsonlite)   # 读 jplace(JSON)
library(ape)        # 读/写 Newick & 处理 phylo
library(treeio)     # treedata/phylo 与 jplace 支持
library(ggtree)     # 系统发育树可视化
library(dplyr)      # 数据操作
library(ggplot2)    # 作图
library(colorspace) # 连续色带

# =========================
# 输入文件路径
# =========================
jplace_file <- "/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/output/placer.aln.snp-sites.jplace"

# =========================
# 读取 jplace & 构造 phylo
# =========================
# 说明：
# 1) 优先使用 treeio::read.jplace（更稳、API 正式）
# 2) 若需要从原始 JSON 列表恢复，则做降级处理（去除 {edge_number} 注释再读 Newick）
# 3) 放置信息使用 treeio 内部提取器作为兜底（实测最健壮）
# =========================

# 优先：正规 API
jobj <- tryCatch(treeio::read.jplace(jplace_file), error = function(e) NULL)
if (!is.null(jobj)) {
  # read.jplace 返回的是 "jplace" S4 对象，含 phylo 与 placements 信息
  phylo <- jobj@phylo

  # placements：尽量用官方导出方法；若不可用，降级到内部提取器
  placements <- tryCatch({
    # 新版 treeio 会把放置信息整理进 as_tibble(jobj, what="placement")
    suppressWarnings(as_tibble(jobj, what = "placement"))
  }, error = function(e) NULL)

  if (is.null(placements) || !"like_weight_ratio" %in% names(placements)) {
    # 兜底：内部提取器（在不同版本 treeio 中最兼容）
    jlist <- jsonlite::fromJSON(jplace_file)
    placements <- getFromNamespace("extract.placement", "treeio")(jlist, phylo)
  }

} else {
  # 降级路径：直接读 JSON + 自行构造 phylo（去掉 {…} 注释）
  jlist <- jsonlite::fromJSON(jplace_file)

  # jplace 树的 Newick 常带 {edge_num} 注解，ape::read.tree 无法直接解析——先去掉
  newick_clean <- gsub("\\{[^}]*\\}", "", jlist$tree)
  phylo <- ape::read.tree(text = newick_clean)

  # 放置信息用 treeio 的内部函数兜底（该函数接受 jplace 列表与 phylo）
  placements <- getFromNamespace("extract.placement", "treeio")(jlist, phylo)
}
placements |> as_tibble() |> write.csv("/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/output/df_pplacer_likelihood.csv")
# 简要健壮性检查
stopifnot(inherits(phylo, "phylo"))
stopifnot(is.data.frame(placements), all(c("name","node","like_weight_ratio") %in% names(placements)))

# =========================
# 过滤：每条序列保留 like_weight_ratio 最大的放置
# =========================
filtered_placements <- placements %>%
  group_by(name) %>%
  filter(like_weight_ratio == max(like_weight_ratio, na.rm = TRUE)) %>%
  ungroup()

# =========================
# 分布对比数据（筛选前 vs 筛选后）
# =========================
dat <- bind_rows(
  transmute(placements, group = "Before",  likelihood_weight_ratio = like_weight_ratio),
  transmute(filtered_placements, group = "Filtered", likelihood_weight_ratio = like_weight_ratio)
)

# =========================
# 直方图：like_weight_ratio 分布
# =========================
p1 <- ggplot(dat, aes(x = likelihood_weight_ratio, fill = group)) +
  geom_histogram(position = "dodge", color = "white", linewidth = 0.5, binwidth = 0.1) +
  scale_fill_manual(values = c("Before" = "#f08c8c", "Filtered" = "#8dc5fe")) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  labs(x = "Likelihood Weight Ratio", y = "Number of Placements") +
  theme_bw()

# =========================
# treedata：把节点放置数合并入树
# =========================
tree_td <- as.treedata(phylo)              # treedata
td_tbl  <- as_tibble(tree_td)              # 节点表（含 node/label/isTip 等）

node_counts <- filtered_placements %>%
  count(node, name = "nplace")

# 合并并回写 treedata 的数据表
tree_td@data <- td_tbl %>%
  left_join(node_counts, by = "node") %>%
  mutate(nplace = ifelse(is.na(nplace), 0L, nplace))

# 供快速查看：按放置数降序
as_tibble(tree_td) %>% arrange(desc(nplace)) %>% head(10)

# =========================
# 圆形树：节点颜色/大小 = 放置数
# =========================
max_n <- max(tree_td@data$nplace %||% 0L, na.rm = TRUE)  # 自动上限
# 避免 0 上限导致 scale 报错
max_n <- if (is.finite(max_n) && max_n > 0) max_n else 1

p2 <- ggtree(
  tree_td,
  layout = "circular",
  branch.length = "none",
  lwd = 1,
  aes(color = nplace, size = nplace)
) +
  colorspace::scale_color_continuous_sequential(palette = "Sunset", limits = c(0, max_n)) +
  scale_size_continuous(range = c(0.1, 2), limits = c(0, max_n),
                        guide = guide_legend(reverse = TRUE)) +
  theme(
    legend.position = "right",
    legend.box = "horizontal",
    legend.box.just = "left"
  )

# 输出图对象（在交互式会话里会自动渲染）
p1
p2
