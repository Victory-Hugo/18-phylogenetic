library(tidyplots)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(colorspace)
library(ragg)
library(scales)
library(ggtext)
library(ggnewscale)

#* =====参数解析=====
args_raw <- commandArgs(trailingOnly = TRUE)
args <- setNames(args_raw[c(FALSE, TRUE)], sub("^--", "", args_raw[c(TRUE, FALSE)]))

frequency_tsv <- args[["frequency-tsv"]]
summary_tsv <- args[["summary-tsv"]]
color_tsv <- args[["color-tsv"]]
label_color_tsv <- args[["label-color-tsv"]]
out_dir <- args[["output-figure"]]
afm_dir <- args[["afm-dir"]]
width_per_population <- as.numeric(args[["width-per-population"]])
panel_height <- as.numeric(args[["panel-height"]])
legend_max_rows <- as.integer(args[["legend-max-rows"]])
dpi <- as.integer(args[["dpi"]])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#* =====字体注册=====
# 注册项目内的 Arial 字体度量，使 grDevices::pdf() 写出的 BaseFont 为 Arial，
# 每个标签都是可在 Illustrator 中单独选中编辑的文本。
# 度量文件缺失时回落到默认字体，不让整条流程因字体而失败。
afm_files <- file.path(afm_dir, c("Arial.afm", "Arial-Bold.afm",
                                  "Arial-Italic.afm", "Arial-BoldItalic.afm"))
plot_font <- if (all(file.exists(afm_files))) {
  grDevices::pdfFonts(Arial = grDevices::Type1Font("Arial", afm_files))
  "Arial"
} else {
  message("未找到 Arial 字体度量，回落到默认字体: ", afm_dir)
  ""
}

#* =====读取数据=====
df1 <- read_tsv(frequency_tsv, show_col_types = FALSE)
df2 <- read_tsv(summary_tsv, show_col_types = FALSE)
print(knitr::kable(head(df1)))

global_level <- df2 |>
  filter(Method == "Global uniform depth", Selected) |>
  pull(`Parameter Value`)
best_lambda <- df2 |>
  filter(Method == "Branch-adaptive depth", Selected) |>
  pull(`Parameter Value`)

scheme_levels <- c("Global uniform depth", "Branch-adaptive depth")
scheme_titles <- c(
  paste0("Global uniform depth: Level ", global_level),
  paste0("Branch-adaptive depth: lambda = ", signif(best_lambda, 3))
)

df_pop <- df1 |>
  distinct(Population, `Population Order`, `Group Label`, `Sample Size`) |>
  arrange(`Population Order`)
pop_levels <- df_pop$Population
n_pop <- length(pop_levels)

#* =====群体标签配色=====
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
         "'>■</span>")
)

# 无群体标签色块时，普通文本不应经 gridtext/Markdown 渲染；后者在旋转
# PDF 文字时会把每个轴标签拆成多个片段。只有实际含 HTML 色块时才使用它。
axis_text_x <- if (all(is.na(df_pop$`Group Label`))) {
  element_text(angle = 90, hjust = 1, vjust = 0.5,
               family = plot_font, size = 7)
} else {
  element_markdown(angle = 90, hjust = 1, vjust = 0.5,
                   family = plot_font, size = 7)
}

#* =====系统发育配色=====
# 与步骤 1 共用主干色表，后代在同色系内由深到浅，祖先残留类别用低饱和浅色
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

hap_palette <- setNames(palette_tbl$Color, palette_tbl$Haplogroup)
hap_levels <- palette_tbl$Haplogroup

#* =====两种方案的堆叠柱状图=====
plot_list <- list()

for (i in seq_along(scheme_levels)) {
  df3 <- df1 |>
    filter(Scheme == scheme_levels[i]) |>
    mutate(
      Population = factor(Population, levels = pop_levels, labels = pop_labels),
      Haplogroup = factor(Haplogroup, levels = hap_levels)
    ) |>
    arrange(Population, Haplogroup)

  n_cat <- nlevels(droplevels(df3$Haplogroup))

  p <- df3 |>
    tidyplot(x = Population, y = Frequency, color = Haplogroup) |>
    add_barstack_relative(width = 0.85, reverse = TRUE) |>
    adjust_colors(new_colors = hap_palette) |>
    adjust_font(family = plot_font, face = "plain", fontsize = 7) |>
    adjust_x_axis(rotate_labels = 90) |>
    adjust_x_axis_title("") |>
    adjust_y_axis(labels = scales::percent) |>
    adjust_y_axis_title(paste0(scheme_titles[i], " (", n_cat, " categories)")) |>
    adjust_legend_title(scheme_levels[i]) |>
    adjust_legend_position("right") |>
    adjust_size(width = n_pop * width_per_population * 25.4,
                height = panel_height * 25.4) |>
    remove_x_axis_ticks()

  p <- p + guides(
    fill = guide_legend(ncol = ceiling(n_cat / legend_max_rows),
                        byrow = FALSE, order = 1,
                        theme = theme(
                          legend.key.size = unit(2.6, "mm"),
                          legend.text = element_text(family = plot_font,
                                                     size = 5),
                          legend.title = element_text(family = plot_font,
                                                      size = 6.5))),
    color = "none"
  ) + theme(
    legend.justification = c(0, 0.5),
    legend.key.spacing.x = unit(1.5, "mm"),
    axis.text.x = axis_text_x
  )

  plot_list[[i]] <- p
}

#* =====群体标签图例=====
# 数据未提供分组标签时（Group Label 全为空），跳过这组图例
if (length(label_levels) > 0) {
  df_lab <- data.frame(
    Population = factor(pop_labels[1], levels = pop_labels),
    Frequency = 0,
    Label = factor(label_levels, levels = label_levels)
  )

  plot_list[[1]] <- plot_list[[1]] +
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
} else {
  message("未检测到群体分组标签，跳过群体标签图例")
}

plot_list[[1]] <- plot_list[[1]] |> remove_x_axis_labels()

p3 <- wrap_plots(plot_list, ncol = 1) +
  plot_annotation(
    title = paste("Y-chromosomal haplogroup composition at the",
                  "cross-validated optimal phylogenetic resolution"),
    theme = theme(
      plot.title = element_text(family = plot_font, face = "plain", size = 12,
                                hjust = 0),
      plot.margin = margin(6, 6, 6, 6, "mm")
    )
  ) &
  theme(text = element_text(family = plot_font, face = "plain"),
        plot.margin = margin(1, 2, 1, 2, "mm"))

#* =====交叉验证诊断图=====
cv_list <- list()

for (i in seq_along(scheme_levels)) {
  df4 <- df2 |>
    filter(Method == scheme_levels[i]) |>
    mutate(Lower = `CV Mean Log Loss` - `CV SE`,
           Upper = `CV Mean Log Loss` + `CV SE`)
  x_title <- if (i == 1) "Phylogenetic level" else "Complexity penalty lambda"

  p <- df4 |>
    tidyplot(x = `Parameter Value`, y = `CV Mean Log Loss`) |>
    add_line(linewidth = 0.4) |>
    add_data_points(size = 1) |>
    adjust_font(family = plot_font, face = "plain", fontsize = 7) |>
    adjust_x_axis_title(x_title) |>
    adjust_y_axis_title("Cross-validated log loss") |>
    adjust_size(width = 70, height = 45) |>
    remove_legend()

  p <- p +
    geom_errorbar(data = df4, inherit.aes = FALSE,
                  aes(x = `Parameter Value`, ymin = Lower, ymax = Upper),
                  width = 0, linewidth = 0.3) +
    geom_point(data = filter(df4, `Best Predictive`), inherit.aes = FALSE,
               aes(x = `Parameter Value`, y = `CV Mean Log Loss`),
               size = 2.2, shape = 21, fill = NA, colour = "#1F77B4",
               stroke = 0.6) +
    geom_point(data = filter(df4, Selected), inherit.aes = FALSE,
               aes(x = `Parameter Value`, y = `CV Mean Log Loss`),
               size = 1.4, colour = "#D62728")
  if (i == 2) p <- p + scale_x_log10()

  cv_list[[i]] <- p
}

df5 <- df2 |>
  filter(Method == "Global uniform depth") |>
  mutate(Lower = `Resolved Rate`, Upper = `Resolved Rate`)

p6 <- df5 |>
  tidyplot(x = `Parameter Value`, y = `Resolved Rate`) |>
  add_line(linewidth = 0.4) |>
  add_data_points(size = 1) |>
  adjust_font(family = plot_font, face = "plain", fontsize = 7) |>
  adjust_x_axis_title("Phylogenetic level") |>
  adjust_y_axis_title("Resolved sample rate") |>
  adjust_size(width = 70, height = 45) |>
  remove_legend()

p7 <- wrap_plots(c(cv_list, list(p6)), nrow = 1) +
  plot_annotation(
    title = paste("Selection of the optimal phylogenetic resolution:",
                  "blue circle marks the lowest cross-validated log loss,",
                  "red dot the one-standard-error choice"),
    theme = theme(
      plot.title = element_text(family = plot_font, face = "plain", size = 10,
                                hjust = 0),
      plot.margin = margin(6, 6, 6, 6, "mm")
    )
  ) &
  theme(text = element_text(family = plot_font, face = "plain"))

#* =====导出=====
n_cat_max <- df1 |> distinct(Scheme, Haplogroup) |> count(Scheme) |> pull(n) |> max()
legend_cols_max <- ceiling(n_cat_max / legend_max_rows)
label_width <- max(nchar(unique(df1$Haplogroup))) * 0.045 + 0.25
combo_w <- n_pop * width_per_population + legend_cols_max * label_width + 1.8
combo_h <- panel_height * length(plot_list) + 4.5

# 用 grDevices::pdf() 而非 cairo_pdf()，避免连续文字被合并进同一个文本对象，
# 保证每个标签在 Illustrator 中可以单独选中编辑
grDevices::pdf(
  file.path(out_dir, "⭐2-4-Optimal-Depth-Frequency.pdf"),
  width = combo_w, height = combo_h, onefile = TRUE,
  family = plot_font, useDingbats = FALSE, encoding = "WinAnsi.enc")
print(p3)
print(p7)
grDevices::dev.off()

ggsave(file.path(out_dir, "⭐2-4-Optimal-Depth-Frequency.png"), p3,
       width = combo_w, height = combo_h, dpi = dpi, device = agg_png,
       limitsize = FALSE)
cat("图形写出目录:", out_dir, "\n")
