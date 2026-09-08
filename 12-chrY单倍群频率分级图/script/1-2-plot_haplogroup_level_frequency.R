library(tidyplots)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(colorspace)
library(ragg)
library(ggtext)
library(ggnewscale)
library(scales)

#* =====参数解析=====
args_raw <- commandArgs(trailingOnly = TRUE)
args <- setNames(args_raw[c(FALSE, TRUE)], sub("^--", "", args_raw[c(TRUE, FALSE)]))

input_tsv <- args[["input-tsv"]]
color_tsv <- args[["color-tsv"]]
label_color_tsv <- args[["label-color-tsv"]]
out_dir <- args[["output-figure"]]
width_per_population <- as.numeric(args[["width-per-population"]])
panel_height <- as.numeric(args[["panel-height"]])
legend_max_rows <- as.integer(args[["legend-max-rows"]])
legend_max_categories <- as.integer(args[["legend-max-categories"]])
dpi <- as.integer(args[["dpi"]])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_font <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else ""

#* =====读取数据=====
df1 <- read_tsv(input_tsv, show_col_types = FALSE)
print(knitr::kable(head(df1)))

# 层级按深度数值排序，避免 Level 10 排到 Level 2 之前
# 类别过多的层级把完整图例移到单独页面，避免图例挤占柱状区
level_stats <- df1 |>
  distinct(Level, `Level Depth`, Haplogroup, `Major Haplogroup`) |>
  group_by(Level, `Level Depth`) |>
  summarise(n_cat = n_distinct(Haplogroup),
            n_major = n_distinct(`Major Haplogroup`), .groups = "drop") |>
  arrange(`Level Depth`) |>
  mutate(external_legend = n_cat > legend_max_categories,
         legend_cols = ceiling(n_cat / legend_max_rows))

inline_levels <- level_stats$Level[!level_stats$external_legend]
max_legend_cols <- max(c(level_stats$legend_cols[!level_stats$external_legend],
                         1))
label_chars <- max(c(nchar(df1$Haplogroup[df1$Level %in% inline_levels]), 4))
legend_width <- max_legend_cols * (label_chars * 1.15 + 6)
df_pop <- df1 |>
  distinct(Population, `Population Order`, `Group Label`, `Sample Size`) |>
  arrange(`Population Order`)
pop_levels <- df_pop$Population
n_pop <- length(pop_levels)

#* =====群体标签配色=====
# 每个分组标签一个方块色，色表可在 conf 中修改；表内没有的标签自动补色
df_label <- read_tsv(label_color_tsv, show_col_types = FALSE)
label_levels <- df_label$Label[df_label$Label %in% df_pop$`Group Label`]
label_extra <- setdiff(unique(na.omit(df_pop$`Group Label`)), label_levels)
label_levels <- c(label_levels, sort(label_extra))
label_color <- setNames(df_label$Color, df_label$Label)[label_levels]
if (length(label_extra) > 0) {
  label_color[label_extra] <- grDevices::colorRampPalette(
    as.character(colors_discrete_rainbow))(length(label_extra))
}

# 群体名后缀一个标签色块，颜色标出该群体所属分组
pop_labels <- ifelse(
  is.na(df_pop$`Group Label`),
  paste0(df_pop$Population, " (N = ", df_pop$`Sample Size`, ")"),
  paste0(df_pop$Population, " (N = ", df_pop$`Sample Size`, ") ",
         "<span style='color:", label_color[df_pop$`Group Label`],
         "'>\u25a0</span>")
)

#* =====系统发育配色=====
# 主干单倍群的基础色由 conf/color.tsv 指定，可自由修改；
# 其后代在同色系内由深到浅区分，祖先残留类别（带星号）用同色系低饱和浅色
df_color <- read_tsv(color_tsv, show_col_types = FALSE)
major_color <- setNames(df_color$Color, df_color$Haplogroup)

missing_major <- setdiff(unique(df1$`Major Haplogroup`), names(major_color))
if (length(missing_major) > 0) {
  stop("以下主干单倍群缺少配色，请补充到 ", color_tsv, ": ",
       paste(missing_major, collapse = ", "))
}

palette_tbl <- df1 |>
  distinct(Haplogroup, `Major Haplogroup`, Resolution, `Category Order`) |>
  arrange(`Category Order`) |>
  group_by(`Major Haplogroup`, Resolution) |>
  mutate(rank_in_group = row_number(), n_in_group = n()) |>
  ungroup() |>
  mutate(
    base = major_color[`Major Haplogroup`],
    shade = ifelse(n_in_group == 1, 0,
                   (rank_in_group - 1) / pmax(n_in_group - 1, 1)),
    Color = ifelse(
      Resolution == "Resolved",
      lighten(base, amount = -0.28 + shade * 0.58),
      desaturate(lighten(base, amount = 0.62 + shade * 0.18), amount = 0.62)
    )
  )

stopifnot(!any(is.na(palette_tbl$Color)))
hap_palette <- setNames(palette_tbl$Color, palette_tbl$Haplogroup)
hap_levels <- palette_tbl$Haplogroup

#* =====逐层级绘图=====
# plot_list：所有样本归入祖先残留类别，柱高恒为 100%
# plot_list2：只画能解析到当前层级的类别，未解析部分在柱顶留空
plot_list <- list()
plot_list2 <- list()

level_names <- level_stats$Level

for (lv in level_names) {
  stats <- level_stats[level_stats$Level == lv, ]
  df2 <- df1 |>
    filter(Level == lv) |>
    mutate(
      Population = factor(Population, levels = pop_levels, labels = pop_labels),
      Haplogroup = factor(Haplogroup, levels = hap_levels)
    ) |>
    arrange(Population, Haplogroup)

  n_cat <- stats$n_cat
  legend_cols <- stats$legend_cols

  p <- df2 |>
    tidyplot(x = Population, y = Frequency, color = Haplogroup) |>
    add_barstack_relative(width = 0.85, reverse = TRUE) |>
    adjust_colors(new_colors = hap_palette) |>
    adjust_font(family = plot_font, face = "plain", fontsize = 7) |>
    adjust_x_axis(rotate_labels = 90) |>
    adjust_x_axis_title("") |>
    adjust_y_axis(labels = scales::percent) |>
    adjust_y_axis_title(paste0(lv, " (", n_cat, " categories)")) |>
    adjust_legend_title(sub("Level", "Level ", lv)) |>
    adjust_legend_position("right") |>
    adjust_size(width = n_pop * width_per_population * 25.4,
                height = panel_height * 25.4) |>
    remove_x_axis_ticks()

  # 各层级图例区等宽并左对齐，保证所有 panel 的柱状区严格对齐
  p <- p + guides(
    fill = guide_legend(ncol = legend_cols, byrow = FALSE, order = 1,
                        keywidth = unit(3, "mm"), keyheight = unit(3, "mm")),
    color = "none"
  ) + theme(
    legend.justification = c(0, 0.5),
    legend.key.spacing.x = unit(1.5, "mm"),
    axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5,
                                   family = plot_font, size = 7),
    legend.text = element_text(family = plot_font, size = 7),
    legend.title = element_text(family = plot_font, size = 7)
  )

  if (stats$external_legend) p <- p |> remove_legend()

  plot_list[[lv]] <- p

  df3 <- df2 |> filter(Resolution == "Resolved")
  n_cat2 <- nlevels(droplevels(df3$Haplogroup))

  p2 <- df3 |>
    tidyplot(x = Population, y = Frequency, color = Haplogroup) |>
    add_barstack_absolute(width = 0.85, reverse = TRUE) |>
    adjust_colors(new_colors = hap_palette) |>
    adjust_font(family = plot_font, face = "plain", fontsize = 7) |>
    adjust_x_axis(rotate_labels = 90) |>
    adjust_x_axis_title("") |>
    adjust_y_axis(limits = c(0, 100), labels = label_percent(scale = 1)) |>
    adjust_y_axis_title(paste0(lv, " (", n_cat2, " categories)")) |>
    adjust_legend_title(sub("Level", "Level ", lv)) |>
    adjust_legend_position("right") |>
    adjust_size(width = n_pop * width_per_population * 25.4,
                height = panel_height * 25.4) |>
    remove_x_axis_ticks()

  p2 <- p2 + guides(
    fill = guide_legend(ncol = ceiling(n_cat2 / legend_max_rows),
                        byrow = FALSE, order = 1,
                        keywidth = unit(3, "mm"), keyheight = unit(3, "mm")),
    color = "none"
  ) + scale_x_discrete(drop = FALSE) + theme(
    legend.justification = c(0, 0.5),
    legend.key.spacing.x = unit(1.5, "mm"),
    axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5,
                                   family = plot_font, size = 7),
    legend.text = element_text(family = plot_font, size = 7),
    legend.title = element_text(family = plot_font, size = 7)
  )

  if (stats$external_legend) p2 <- p2 |> remove_legend()

  plot_list2[[lv]] <- p2
}

#* =====群体标签图例=====
# 借用首个面板右侧的图例空白，追加一组群体标签色块说明
# 数据未提供分组标签时（Group Label 全为空），跳过这组图例
if (length(label_levels) > 0) {
  df_lab <- data.frame(
    Population = factor(pop_labels[1], levels = pop_labels),
    Frequency = 0,
    Label = factor(label_levels, levels = label_levels)
  )

  for (i in c(1, 2)) {
    target <- if (i == 1) plot_list else plot_list2
    target[[1]] <- target[[1]] +
      new_scale_fill() +
      geom_tile(data = df_lab, aes(x = Population, y = Frequency, fill = Label),
                width = 0, height = 0, inherit.aes = FALSE) +
      scale_fill_manual(values = label_color, name = "Population group") +
      guides(fill = guide_legend(
        ncol = 1, order = 2,
        theme = theme(legend.key.size = unit(3, "mm"),
                      legend.text = element_text(family = plot_font, size = 5.5),
                      legend.title = element_text(family = plot_font,
                                                  size = 6.5))))
    if (i == 1) plot_list <- target else plot_list2 <- target
  }
} else {
  message("未检测到群体分组标签，跳过群体标签图例")
}

#* =====纵向组合图=====
# 上方各 panel 隐藏重复的群体名称，仅最下方保留完整横轴
for (i in seq_along(plot_list)[-length(plot_list)]) {
  plot_list[[i]] <- plot_list[[i]] |> remove_x_axis_labels()
  plot_list2[[i]] <- plot_list2[[i]] |> remove_x_axis_labels()
}

p3 <- wrap_plots(plot_list, ncol = 1) +
  plot_annotation(
    title = "Y-chromosomal haplogroup composition across phylogenetic levels",
    theme = theme(
      plot.title = element_text(family = plot_font, face = "plain", size = 12,
                                hjust = 0),
      plot.margin = margin(6, 6, 6, 6, "mm")
    )
  ) &
  theme(text = element_text(family = plot_font, face = "plain"),
        plot.margin = margin(1, 2, 1, 2, "mm"))

combo_w <- n_pop * width_per_population + legend_width / 25.4 + 1.6
combo_h <- panel_height * length(plot_list) + 4.5

p4 <- wrap_plots(plot_list2, ncol = 1) +
  plot_annotation(
    title = paste("Y-chromosomal haplogroup composition across phylogenetic",
                  "levels, unresolved haplogroups left blank"),
    theme = theme(
      plot.title = element_text(family = plot_font, face = "plain", size = 12,
                                hjust = 0),
      plot.margin = margin(6, 6, 6, 6, "mm")
    )
  ) &
  theme(text = element_text(family = plot_font, face = "plain"),
        plot.margin = margin(1, 2, 1, 2, "mm"))

#* =====外置图例页=====
# 类别过多的层级，完整图例单独成页，颜色与柱状图完全一致
legend_pages <- list()
rows_per_column <- floor((combo_h - 1.5) * 25.4 / 3.4)

for (lv in level_stats$Level[level_stats$external_legend]) {
  cats <- hap_levels[hap_levels %in% df1$Haplogroup[df1$Level == lv]]
  df6 <- data.frame(x = 1, y = 1,
                    Haplogroup = factor(cats, levels = cats))

  legend_pages[[lv]] <- ggplot(df6, aes(x, y, fill = Haplogroup)) +
    geom_tile(width = 0, height = 0) +
    scale_fill_manual(values = hap_palette[cats],
                      name = paste0(lv, " (", length(cats), " categories)")) +
    guides(fill = guide_legend(
      ncol = ceiling(length(cats) / rows_per_column), byrow = FALSE,
      theme = theme(legend.key.size = unit(2.6, "mm"),
                    legend.text = element_text(family = plot_font, size = 5),
                    legend.title = element_text(family = plot_font,
                                                size = 9)))) +
    theme_void() +
    theme(legend.position = "left",
          legend.justification = c(0, 1),
          legend.key.spacing.x = unit(1.5, "mm"),
          plot.margin = margin(6, 6, 6, 6, "mm"),
          text = element_text(family = plot_font, face = "plain"))
}

# 两个版本与外置图例合并为同一 PDF 的多页
grDevices::cairo_pdf(
  file.path(out_dir, "⭐1-2-Haplogroup-Frequency-All-Levels.pdf"),
  width = combo_w, height = combo_h, onefile = TRUE)
print(p3)
print(p4)
for (lp in legend_pages) print(lp)
grDevices::dev.off()

png_dpi <- min(dpi, floor(50000 / max(combo_w, combo_h)),
               floor(sqrt(2e8 / (combo_w * combo_h))))
if (png_dpi < dpi) {
  message("画布 ", round(combo_w, 1), " x ", round(combo_h, 1),
          " inch 超过位图上限，PNG 分辨率下调为 ", png_dpi, " dpi")
}

ggsave(file.path(out_dir, "⭐1-2-Haplogroup-Frequency-All-Levels.png"), p3,
       width = combo_w, height = combo_h, dpi = png_dpi, device = agg_png,
       limitsize = FALSE)
ggsave(file.path(out_dir, "1-2-Haplogroup-Frequency-Unresolved-Blank.png"), p4,
       width = combo_w, height = combo_h, dpi = png_dpi, device = agg_png,
       limitsize = FALSE)
if (length(legend_pages) > 0) {
  ggsave(file.path(out_dir, "1-2-Haplogroup-Frequency-Legends.png"),
         wrap_plots(legend_pages, nrow = 1),
         width = combo_w, height = combo_h, dpi = png_dpi, device = agg_png,
         limitsize = FALSE)
  message("类别数超过 ", legend_max_categories, " 的层级共 ",
          length(legend_pages), " 个，完整图例已单独成页")
}

write_tsv(df1, file.path(out_dir, "⭐1-2-Haplogroup-Frequency-All-Levels.tsv"))

cat("图形写出目录:", out_dir, "\n")
