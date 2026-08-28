#* =====依赖与参数=====
library(readr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(grid)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
cluster_file <- args[1] #* ⭐2-Cluster-And-Similarity.tsv
fig_dir <- args[2]
palette_spec <- if (length(args) >= 3) args[3] else "diverging"

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "palette.R"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

#* =====读取数据=====
# 两张子图：水平分数（平均 log Ne 与波动幅度）与变化分数（净变化、速率、扩张与收缩强度）。
# 合成单一分数行不通——实测那一个主成分的载荷几乎全落在变化类特征上，图里根本没有水平
# 信息。颜色用逐时间段标准化后的分数：同一列内比较"谁偏高、谁偏低"，不跨列读绝对高低。
dat1 <- read_tsv(cluster_file, show_col_types = FALSE) %>%
  filter(Quantity == "Standardized Composite Score") %>%
  mutate(Value = as.numeric(Value),
         `Window Start (Years)` = as.numeric(`Window Start (Years)`)) %>%
  select(Category, Group, `Score Component`, `Time Window`,
         `Window Start (Years)`, Score = Value)
if (!nrow(dat1)) stop("结果表里没有标准化分数行", call. = FALSE)

order_win <- dat1 %>%
  distinct(`Time Window`, `Window Start (Years)`) %>%
  arrange(`Window Start (Years)`) %>%
  pull(`Time Window`)
components <- sort(unique(dat1$`Score Component`))

to_matrix <- function(d) {
  m <- d %>%
    select(Category, `Time Window`, Score) %>%
    pivot_wider(names_from = `Time Window`, values_from = Score) %>%
    as.data.frame()
  rownames(m) <- m$Category
  m$Category <- NULL
  as.matrix(m)[, intersect(order_win, colnames(m)), drop = FALSE]
}
mats <- lapply(components, function(cp) to_matrix(dat1 %>% filter(`Score Component` == cp)))
names(mats) <- components

#* =====行序=====
# 行聚类只在"所有类别都可靠"的时间列上做：缺格多时按两两可用列算相关，距离由少数几列
# 决定，聚出来的结构不稳也不可解释
joint <- do.call(rbind, mats)
common_cols <- colnames(mats[[1]])[colSums(is.na(joint)) == 0]
row_tree <- NULL
if (length(common_cols) >= 3 && nrow(mats[[1]]) > 2) {
  block <- do.call(cbind, lapply(mats, function(m) m[, common_cols, drop = FALSE]))
  row_tree <- hclust(dist(block, method = "euclidean"), method = "average")
} else {
  message("共同可靠的时间列不足 3 个，行聚类已关闭")
}

anno_row <- dat1 %>% distinct(Category, Group) %>% as.data.frame()
rownames(anno_row) <- anno_row$Category
anno_row$Category <- NULL
group_levels <- sort(unique(anno_row$Group))
pal_group <- colorRampPalette(pal_discrete)(length(group_levels))
names(pal_group) <- group_levels

# 色标以 0 为中心并取 98 分位为上下界：个别极端格子会把整张图压成一片同色
lim <- as.numeric(quantile(abs(joint), 0.98, na.rm = TRUE))
if (!is.finite(lim) || lim <= 0) lim <- max(abs(joint), na.rm = TRUE)
brk <- seq(-lim, lim, length.out = 101)
heat_cols <- resolve_palette(palette_spec)

step <- max(1, ceiling(ncol(mats[[1]]) / 20))
lab_col <- ifelse(seq_len(ncol(mats[[1]])) %% step == 1, colnames(mats[[1]]), "")

grobs <- lapply(components, function(cp) {
  m <- pmin(pmax(mats[[cp]], -lim), lim)
  ph <- pheatmap(m,
    #* 两个子图共用同一棵树：行序一致才谈得上上下对照
    cluster_rows = if (is.null(row_tree)) FALSE else row_tree,
    treeheight_row = 24,
    cluster_cols = FALSE,          #* 列是时间，顺序不可打乱
    cellwidth = 10, cellheight = 10,   #* 正方形格子
    angle_col = 90,
    fontsize = 8,
    fontsize_row = 8,
    fontsize_col = 6,
    labels_col = lab_col,
    annotation_row = anno_row,
    annotation_colors = list(Group = pal_group),
    annotation_legend = FALSE,   #* 图例统一放到整版图右侧
    legend = FALSE,
    main = cp,
    breaks = brk,
    color = heat_cols,
    na_col = na_colour,
    border_color = "#FFFFFF",   # 纯白色格子描边
    silent = TRUE
  )
  plain_pheatmap(ph)$gtable
})

stem <- file.path(fig_dir, "2-Composite-Score-Heatmap")
save_panels(grobs, stem, n_col = 1, panel_width = 7.6, panel_height = 2.5,
            title = grid::textGrob("Composite trajectory scores across the whole time axis",
                                   gp = grid::gpar(fontfamily = plot_font_family,
                                                   fontface = 1, fontsize = 12)),
            legend = shared_legend(heat_cols, lim, "Standardized score",
                                   group_levels, pal_group),
            legend_width = 1.8)
write_tsv(dat1 %>% arrange(`Score Component`, Category, `Window Start (Years)`),
          paste0(stem, ".tsv"))
