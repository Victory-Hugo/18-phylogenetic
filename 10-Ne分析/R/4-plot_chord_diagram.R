#* =====依赖与参数=====
library(readr)
library(dplyr)
library(tibble)
library(igraph)
library(ggraph)
library(ggforce)
library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
cluster_file <- args[1] #* ⭐2-Cluster-And-Similarity.tsv
fig_dir <- args[2]
threshold <- if (length(args) >= 3) as.numeric(args[3]) else 0.5
n_col <- if (length(args) >= 4) as.integer(args[4]) else 3

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(script_dir, "palette.R"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

#* =====读取数据=====
# 分层弦图：结点按分组排列、外圈色带标出分组，弦的粗细为两个类别的同簇频率
dat1 <- read_tsv(cluster_file, show_col_types = FALSE)
link_quantity <- if (any(dat1$Quantity == "Co-clustering Frequency")) {
  "Co-clustering Frequency"
} else {
  "Similarity"
}
between_label <- "Between groups"

periods <- dat1 %>%
  filter(Quantity == "Node Size") %>%
  mutate(start = as.numeric(`Window Start (Years)`)) %>%
  distinct(`Time Window`, start) %>%
  arrange(start) %>%
  pull(`Time Window`)

# 连边配色在所有子图之间统一，图例只画一份放在整版右侧
all_groups <- dat1 %>%
  filter(Quantity == "Node Size") %>%
  pull(Group) %>%
  unique() %>%
  sort()
group_cols_all <- colorRampPalette(pal_discrete)(length(all_groups))
names(group_cols_all) <- all_groups
link_cols_all <- c(group_cols_all, setNames(na_colour, between_label))

grobs <- list()
all_links <- list()
for (period in periods) {
  nodes <- dat1 %>%
    filter(Quantity == "Node Size", `Time Window` == period) %>%
    transmute(name = Category, group = Group, size = as.numeric(Value)) %>%
    arrange(group, name)
  if (nrow(nodes) < 3) next

  # 弦的颜色：同组用该组的颜色，跨组统一用灰色——否则"用起点的组"会让同一条弦的
  # 归属取决于排序，读者无从解释
  links <- dat1 %>%
    filter(Quantity == link_quantity, `Time Window` == period,
           as.numeric(Value) >= threshold) %>%
    transmute(from = Category, to = `Second Category`, value = as.numeric(Value)) %>%
    left_join(nodes %>% select(name, group_from = group), by = c("from" = "name")) %>%
    left_join(nodes %>% select(name, group_to = group), by = c("to" = "name")) %>%
    mutate(link = ifelse(group_from == group_to, group_from, between_label)) %>%
    # 显式定类型：某个时间段一条边都没有时，ifelse 会给出逻辑型空列，
    # 与其他时间段的字符列拼不到一起
    transmute(from = as.character(from), to = as.character(to),
              value = as.numeric(value), link = as.character(link))

  g <- graph_from_data_frame(links, vertices = nodes, directed = FALSE)
  lay <- create_layout(g, layout = "linear", circular = TRUE)

  group_cols <- group_cols_all
  link_cols <- link_cols_all

  #* 外圈色带的角度直接从结点坐标反推：ggforce 的 arc 角度是"自 12 点起顺时针"，
  #* 与 atan2 给出的数学角方向相反，用下标硬算会让色带整体转过 90°
  n_total <- nrow(lay)
  half_step <- pi / n_total
  gap <- 4 * pi / 180
  ang <- (pi / 2 - atan2(lay$y, lay$x)) %% (2 * pi)
  band <- tibble(group = lay$group, angle = ang, idx = seq_len(n_total)) %>%
    group_by(group) %>%
    summarise(start = min(angle), end = max(angle), n = n(), .groups = "drop") %>%
    # 跨过 12 点的那一组两端角度会分居 0 与 2π，按跨度判断并整体偏移
    mutate(cross = (end - start) > pi & n > 1,
           start2 = ifelse(cross, end, start) - half_step + gap / 2,
           end2 = ifelse(cross, start + 2 * pi, end) + half_step - gap / 2)

  r_nodes <- mean(sqrt(lay$x^2 + lay$y^2))
  r_band0 <- r_nodes + 0.16
  r_band1 <- r_band0 + 0.07
  r_label <- r_band1 + 0.10

  #* 圆周标签：左半圆的文字要转 180° 并右对齐，否则一半标签是倒着的
  lab_df <- lay %>%
    as_tibble() %>%
    mutate(name = nodes$name, group = nodes$group,
           scale = r_label / sqrt(x^2 + y^2),
           x_lab = x * scale, y_lab = y * scale,
           deg = atan2(y, x) * 180 / pi,
           flip = x < 0,
           angle = ifelse(flip, deg + 180, deg),
           hjust = ifelse(flip, 1, 0))

  # 某个时间段可能一条连边都不过阈值，此时不能加边图层：ggraph 会在空数据上报错
  p <- ggraph(lay)
  if (nrow(links)) {
    p <- p +
      geom_edge_arc2(aes(width = value, colour = link), alpha = 0.6, lineend = "round") +
      scale_edge_width(range = c(0.2, 2.6), guide = "none") +
      scale_edge_colour_manual(values = link_cols, name = "Link",
                               breaks = names(link_cols))
  }
  p <- p +
    geom_arc_bar(data = band,
                 aes(x0 = 0, y0 = 0, r0 = r_band0, r = r_band1,
                     start = start2, end = end2, fill = group),
                 colour = NA, alpha = 0.55, inherit.aes = FALSE) +
    scale_fill_manual(values = group_cols, guide = "none") +
    geom_node_point(aes(x = x, y = y, colour = group, size = size), show.legend = FALSE) +
    scale_colour_manual(values = group_cols) +
    scale_size(range = c(3, 8)) +
    geom_text(data = lab_df,
              aes(x = x_lab, y = y_lab, label = name, angle = angle, hjust = hjust),
              colour = "grey20", vjust = 0.5, size = 1.7, show.legend = FALSE) +
    coord_equal(clip = "off", xlim = c(-2.7, 2.7), ylim = c(-2.7, 2.7)) +
    labs(title = period) +
    theme_void(base_family = plot_font_family) +
    guides(colour = "none") +
    theme(text = element_text(face = "plain"),
          plot.title = element_text(face = "plain", hjust = 0.5, size = 8),
          legend.position = "none",
          plot.margin = margin(4, 4, 4, 4),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA))

  grobs[[length(grobs) + 1]] <- ggplot2::ggplotGrob(p)
  all_links[[length(all_links) + 1]] <- links %>% mutate(`Time Window` = period)
}

# 全部时间段排进同一张图：逐段单独出 PDF 时，读者要在几十个文件之间来回翻
stem <- file.path(fig_dir, "4-Chord-Diagram")
save_panels(grobs, stem, n_col = n_col, panel_width = 3.0, panel_height = 3.0,
            title = grid::textGrob("Similarity of trajectories between categories",
                                   gp = grid::gpar(fontfamily = plot_font_family,
                                                   fontface = 1, fontsize = 12)),
            legend = shared_legend(groups = names(link_cols_all),
                                   group_cols = link_cols_all, group_title = "Link"),
            legend_width = 1.6)
write_tsv(bind_rows(all_links), paste0(stem, ".tsv"))
