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
n_col <- if (length(args) >= 4) as.integer(args[4]) else 3

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "palette.R"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

#* =====读取数据=====
# 每个时间段一对子图：左边同簇频率（后验重抽里两个类别落进同一簇的比例），
# 右边特征相似度（1 − 归一化欧氏距离）。同一时段左右并排，"后验上很稳定"与
# "特征上很相似"是否一致一眼可见
QUANTITIES <- c("Co-clustering Frequency", "Similarity")
QUANTITY_TITLE <- c("Co-clustering Frequency" = "Co-clustering frequency",
                    "Similarity" = "Feature similarity")

dat1 <- read_tsv(cluster_file, show_col_types = FALSE) %>%
  filter(Quantity %in% QUANTITIES, `Window Type` == "Time Period") %>%
  mutate(Value = as.numeric(Value),
         start = as.numeric(`Window Start (Years)`)) %>%
  select(Quantity, `Time Window`, start, Category, `Second Category`, Value)
if (!nrow(dat1)) stop("结果表里没有配对量（同簇频率 / 相似度）行", call. = FALSE)

periods <- dat1 %>%
  distinct(`Time Window`, start) %>%
  arrange(start) %>%
  pull(`Time Window`)

# 行列顺序在所有子图之间**只定一次**：按分组、组内按名称。逐子图各聚各的时，同一个
# 格子在不同时间段代表不同的类别对，几十张子图之间根本没法横比
dat2 <- read_tsv(cluster_file, show_col_types = FALSE) %>%
  filter(Quantity == "Cluster Membership", `Window Type` == "Time Period") %>%
  distinct(Category, Group) %>%
  arrange(Group, Category)
cat_order <- dat2$Category

anno_row <- as.data.frame(dat2[, "Group", drop = FALSE])
rownames(anno_row) <- dat2$Category
group_levels <- sort(unique(dat2$Group))
pal_group <- colorRampPalette(pal_discrete)(length(group_levels))
names(pal_group) <- group_levels

# 两个量都落在 0–1，共用一套色标：不设 breaks 时 pheatmap 会逐图自动定标，
# 同一个颜色在不同子图代表不同数值
heat_cols <- resolve_palette(palette_spec)
brk <- seq(0, 1, length.out = length(heat_cols) + 1)

#* =====逐时间段作图=====
grobs <- list()
all_rows <- list()
for (period in periods) {
  pair_grobs <- list()
  for (qty in QUANTITIES) {
    dat3 <- dat1 %>% filter(Quantity == qty, `Time Window` == period)
    if (!nrow(dat3)) next

    # 补齐成完整方阵：只写了上三角，另一半按对称补，对角线按定义置 1；
    # 该时段没进聚类的类别整行整列留 NA，在图上是灰色，不会伪装成"相似度为 0"
    mat1 <- matrix(NA_real_, nrow = length(cat_order), ncol = length(cat_order),
                   dimnames = list(cat_order, cat_order))
    for (i in seq_len(nrow(dat3))) {
      a <- dat3$Category[i]
      b <- dat3$`Second Category`[i]
      mat1[a, b] <- dat3$Value[i]
      mat1[b, a] <- dat3$Value[i]
    }
    present <- rownames(mat1)[rowSums(!is.na(mat1)) > 0]
    diag(mat1)[cat_order %in% present] <- 1

    ph <- pheatmap(
      mat1,
      scale = "none",
      cluster_rows = FALSE,        #* 行列顺序全局固定，不允许逐图重排
      cluster_cols = FALSE,
      cellwidth = 13,
      cellheight = 13,
      angle_col = 90,
      fontsize = 7,
      fontsize_row = 7,
      fontsize_col = 7,
      breaks = brk,
      annotation_row = anno_row,
      annotation_colors = list(Group = pal_group),
      annotation_legend = FALSE,
      legend = FALSE,              #* 图例统一放到整版图右侧
      main = sprintf("%s — %s", period, QUANTITY_TITLE[[qty]]),
      color = heat_cols,
      na_col = na_colour,
      border_color = "#FFFFFF",   # 纯白色格子描边
      silent = TRUE
    )
    pair_grobs[[length(pair_grobs) + 1]] <- plain_pheatmap(ph)$gtable
    all_rows[[length(all_rows) + 1]] <- dat3 %>% mutate(Quantity = qty)
  }
  if (length(pair_grobs) < 2) next
  # 一对矩阵合成一个单元后再排版：拆成两个独立子图的话，一对会被行末截断成两页
  grobs[[length(grobs) + 1]] <- gridExtra::arrangeGrob(grobs = pair_grobs, ncol = 2)
}
if (!length(grobs)) stop("没有任何时间段同时给出两个配对量", call. = FALSE)

#* =====出图=====
# n_col 的语义是"每行几个矩阵子图"，一对占两个位置，因此每行放 n_col %/% 2 对
pairs_per_row <- max(1, n_col %/% 2)
stem <- file.path(fig_dir, "5-Pairwise-Matrix-Heatmap")
save_panels(grobs, stem, n_col = pairs_per_row, panel_width = 8.0, panel_height = 3.4,
            title = grid::textGrob("Pairwise similarity between categories by time period",
                                   gp = grid::gpar(fontfamily = plot_font_family,
                                                   fontface = 1, fontsize = 12)),
            legend = shared_legend(heat_cols, c(0, 1), "Pairwise value",
                                   group_levels, pal_group),
            legend_width = 1.8)
write_tsv(bind_rows(all_rows) %>%
            select(Quantity, `Time Window`, `Window Start (Years)` = start,
                   Category, `Second Category`, Value) %>%
            arrange(Quantity, `Window Start (Years)`, Category, `Second Category`),
          paste0(stem, ".tsv"))
