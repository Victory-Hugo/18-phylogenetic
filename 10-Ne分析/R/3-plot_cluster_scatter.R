#* =====依赖与参数=====
library(readr)
library(dplyr)
library(tidyr)
library(tidyplots)
library(ggrepel)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
cluster_file <- args[1] #* ⭐2-Cluster-And-Similarity.tsv
fig_dir <- args[2]
scatter_method <- if (length(args) >= 3) args[3] else "pca"
show_ellipse <- if (length(args) >= 4) tolower(args[4]) %in% c("true", "1", "yes") else TRUE
n_col <- if (length(args) >= 5) as.integer(args[5]) else 3

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "palette.R"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

#* 自己算椭圆，不用 stat_ellipse：后者要求 n - 1 >= 3，也就是至少 4 个点，3 个点的簇
#* 会被直接跳过。半径取卡方分位（已知协方差下的正态椭圆），而不是小样本的 F 校正——
#* n = 3 时 F 校正给出的 95% 椭圆半径是 6 个标准差，画出来会把整张子图撑爆。
ellipse_paths <- function(d, level = 0.95, segments = 72) {
  rad <- sqrt(stats::qchisq(level, df = 2))
  th <- seq(0, 2 * pi, length.out = segments)
  circle <- cbind(cos(th), sin(th))
  out <- list()
  for (cl in unique(d$Cluster)) {
    sub <- as.matrix(d[d$Cluster == cl, c("Component 1", "Component 2")])
    if (nrow(sub) < 3) next
    ev <- eigen(stats::cov(sub), symmetric = TRUE)
    if (any(ev$values <= 1e-10)) next        #* 点共线时协方差退化，椭圆没有意义
    pts <- rad * circle %*% diag(sqrt(ev$values)) %*% t(ev$vectors)
    out[[cl]] <- data.frame(
      Cluster = cl,
      `Component 1` = pts[, 1] + mean(sub[, 1]),
      `Component 2` = pts[, 2] + mean(sub[, 2]),
      check.names = FALSE)
  }
  if (!length(out)) NULL else do.call(rbind, out)
}

#* =====读取数据=====
dat1 <- read_tsv(cluster_file, show_col_types = FALSE)
axis_label <- function(ratio, i) {
  base <- if (tolower(scatter_method) == "umap") "UMAP dimension" else "PC"
  if (is.na(ratio)) sprintf("%s %d", base, i) else sprintf("%s %d (%.0f%%)", base, i, 100 * ratio)
}

periods <- dat1 %>%
  filter(Quantity == "Scatter Coordinate") %>%
  mutate(start = as.numeric(`Window Start (Years)`)) %>%
  distinct(`Time Window`, start) %>%
  arrange(start) %>%
  pull(`Time Window`)

# 簇配色在所有子图之间统一：簇标签已跨窗对齐过，各子图各调一套色的话，同一个颜色
# 在不同时间段会代表不同的簇
cluster_levels <- dat1 %>%
  filter(Quantity == "Scatter Coordinate") %>%
  pull(Cluster) %>%
  unique() %>%
  sort()
pal_cluster <- colorRampPalette(pal_discrete)(length(cluster_levels))
names(pal_cluster) <- cluster_levels

grobs <- list()
all_rows <- list()
for (period in periods) {
  dat2 <- dat1 %>%
    filter(Quantity == "Scatter Coordinate", `Time Window` == period) %>%
    select(Category, Group, Cluster, Component, Value) %>%
    pivot_wider(names_from = Component, values_from = Value)

  ratio <- dat1 %>%
    filter(Quantity == "Explained Variance", `Time Window` == period) %>%
    arrange(Component) %>%
    pull(Value)
  if (length(ratio) < 2) ratio <- c(NA_real_, NA_real_)

  p <- dat2 %>%
    tidyplot(x = `Component 1`, y = `Component 2`, color = Cluster) |>
    add_data_points(size = 1.6)
  #* 标注每个点对应的 Category：白底半透明，描边与点同色（沿用 tidyplot 的 colour 映射）
  p <- p + ggrepel::geom_label_repel(
    ggplot2::aes(label = Category),
    fill = ggplot2::alpha("white", 0.5),
    label.size = 0.2, label.r = grid::unit(0.5, "mm"),
    label.padding = grid::unit(0.6, "mm"),
    family = plot_font_family, size = 1.5,
    segment.size = 0.2, segment.alpha = 0.6,
    min.segment.length = 0, box.padding = 0.25, point.padding = 0.15,
    max.overlaps = Inf, seed = 1, show.legend = FALSE)

  # 每个成员数不少于 3 的簇都画一个椭圆
  ell <- if (show_ellipse) ellipse_paths(dat2) else NULL
  if (!is.null(ell)) {
    p <- p + ggplot2::geom_path(data = ell, linewidth = 0.3, alpha = 0.55,
                                show.legend = FALSE)
  }
  p <- p |>
    adjust_colors(new_color_scheme(pal_cluster[sort(unique(dat2$Cluster))])) |>
    adjust_font(family = plot_font_family, face = "plain", fontsize = 7) |>
    adjust_x_axis_title(axis_label(ratio[1], 1)) |>
    adjust_y_axis_title(axis_label(ratio[2], 2)) |>
    adjust_title(period) |>
    remove_legend() |>
    adjust_size(width = 58, height = 58)

  grobs[[length(grobs) + 1]] <- ggplot2::ggplotGrob(p)
  all_rows[[length(all_rows) + 1]] <- dat2 %>% mutate(`Time Window` = period)
}

stem <- file.path(fig_dir, "3-Cluster-Scatter")
save_panels(grobs, stem, n_col = n_col, panel_width = 3.5, panel_height = 3.5,
            title = grid::textGrob("Category ordination of trajectory features by time period",
                                   gp = grid::gpar(fontfamily = plot_font_family,
                                                   fontface = 1, fontsize = 12)),
            legend = shared_legend(groups = cluster_levels, group_cols = pal_cluster,
                                   group_title = "Cluster"),
            legend_width = 1.3)
write_tsv(bind_rows(all_rows), paste0(stem, ".tsv"))
