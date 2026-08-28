#* =====依赖与参数=====
library(readr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(grid)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
feature_file <- args[1] #* ⭐1-Feature-Summary.tsv
fig_dir <- args[2]
palette_spec <- if (length(args) >= 3) args[3] else "diverging"
n_col <- if (length(args) >= 4) as.integer(args[4]) else 3

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "palette.R"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

#* =====读取数据=====
# 每个时间段一个子图（类别 × 8 特征），全部子图排进同一张图，每行 n_col 个
dat1 <- read_tsv(feature_file, show_col_types = FALSE) %>%
  filter(`Window Type` == "Time Period", Usable == "Yes")

periods <- dat1 %>%
  distinct(`Time Window`, `Window Start (Years)`) %>%
  arrange(`Window Start (Years)`) %>%
  pull(`Time Window`)
if (!length(periods)) stop("特征表里没有可用的时间段行", call. = FALSE)

heat_cols <- resolve_palette(palette_spec)

# 全部子图共用一套色标：不设 breaks 时 pheatmap 会逐图自动定标，同一个颜色在不同
# 时间段代表不同数值，子图之间根本没法横比
z_lim <- 2.5
brk <- seq(-z_lim, z_lim, length.out = length(heat_cols) + 1)

group_levels <- sort(unique(dat1$Group))
pal_group <- colorRampPalette(pal_discrete)(length(group_levels))
names(pal_group) <- group_levels

# 特征列的顺序**只聚一次**，所有子图共用：逐子图各聚各的时，同一列在不同时间段代表不同
# 特征，几十张子图之间根本没法横比，标签也只能每张重画一遍
pooled <- dat1 %>%
  mutate(cell = paste(Category, `Time Window`)) %>%
  select(cell, Feature, Median) %>%
  pivot_wider(names_from = Feature, values_from = Median) %>%
  as.data.frame()
rownames(pooled) <- pooled$cell
pooled$cell <- NULL
pooled <- scale(as.matrix(pooled))
col_order <- colnames(pooled)[hclust(dist(t(pooled)), method = "average")$order]

grobs <- list()
kept <- character(0)

for (period in periods) {
  dat2 <- dat1 %>%
    filter(`Time Window` == period) %>%
    select(Category, Group, Feature, Median)

  mat1 <- dat2 %>%
    select(Category, Feature, Median) %>%
    pivot_wider(names_from = Feature, values_from = Median) %>%
    as.data.frame()
  rownames(mat1) <- mat1$Category
  mat1$Category <- NULL
  mat1 <- as.matrix(mat1)
  if (nrow(mat1) < 2) next

  # 逐特征 z-score 自己算：某个特征在所有类别上取值相同时标准差为 0，
  # 交给 scale = "column" 会得到整列 NaN 并让 hclust 直接报错
  sd_col <- apply(mat1, 2, sd)
  sd_col[sd_col < 1e-12] <- 1
  mat1 <- sweep(sweep(mat1, 2, colMeans(mat1), "-"), 2, sd_col, "/")

  # 逐特征 z-score 后仍可能有个别极端值，压到色标上下界，否则色标要为它整体拉宽
  mat1 <- pmin(pmax(mat1, -z_lim), z_lim)[, col_order, drop = FALSE]
  colnames(mat1) <- short_feature(colnames(mat1))

  anno_row <- dat2 %>% distinct(Category, Group) %>% as.data.frame()
  rownames(anno_row) <- anno_row$Category
  anno_row$Category <- NULL

  ph <- pheatmap(
      mat1,
      scale = "none",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      cellwidth = 13,
      cellheight = 13,
      angle_col = 90,
      fontsize = 7,
      fontsize_row = 7,
      fontsize_col = 7,
      treeheight_row = 20,
      breaks = brk,
      annotation_row = anno_row,
      annotation_colors = list(Group = pal_group),
      annotation_legend = FALSE,
      legend = FALSE,
      main = period,
      color = heat_cols,
      na_col = na_colour,
      border_color = "#FFFFFF",   # 纯白色格子描边
      silent = TRUE
  )
  grobs[[length(grobs) + 1]] <- plain_pheatmap(ph)$gtable
  kept <- c(kept, period)
}

stem <- file.path(fig_dir, "1-Period-Feature-Heatmap")
save_panels(grobs, stem, n_col = n_col, panel_width = 4.4, panel_height = 2.5,
            title = grid::textGrob("Trajectory features by category and time period",
                                   gp = grid::gpar(fontfamily = plot_font_family,
                                                   fontface = 1, fontsize = 12)),
            legend = shared_legend(heat_cols, z_lim, "Standardized feature",
                                   group_levels, pal_group),
            legend_width = 1.8)
write_tsv(dat1 %>% filter(`Time Window` %in% kept) %>%
            select(Category, Group, `Time Window`, `Window Start (Years)`,
                   Feature, Median),
          paste0(stem, ".tsv"))
