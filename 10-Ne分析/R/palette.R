#* =====项目统一调色盘与图形约定=====
# 所有出图脚本 source() 本文件取色与字体，不使用任何绘图库的默认色板。

pal_discrete_1 <- c("#0072b2", "#56b4e9", "#009e73", "#f0e442", "#e69f00", "#d55e00")
pal_discrete_2 <- c("#8ecae6", "#219ebc", "#023047", "#ffb703", "#fb8500")
pal_discrete <- c(pal_discrete_1, pal_discrete_2)

pal_sequential <- c(
  "#0b0405", "#30203e", "#3e4d93", "#366b9f", "#3488a6", "#36a4ab",
  "#49c1ad", "#60ceac", "#84d8b0", "#c4e9cf", "#def5e5"
)

pal_diverging <- c(
  "#5b53a4", "#456fb1", "#368cbb", "#4fa8af", "#69c3a5", "#8ad0a4",
  "#addda3", "#c9e99d", "#e6f598", "#fefdbc", "#fef0a6", "#fdc877",
  "#f88f52", "#e55848", "#bc2349", "#a20643"
)

plot_font_family <- "Arial"
na_colour <- "grey88"          #* 未通过覆盖度判据的格子一律留灰

#* 图例不带刻度线；ggplot2 4.x 才有 legend.ticks，低版本静默跳过
drop_legend_ticks <- function(p) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    p + ggplot2::theme(legend.ticks = ggplot2::element_blank())
  } else {
    p
  }
}

#* pheatmap 的所有文字压回 regular 字重
plain_pheatmap <- function(ph) {
  for (i in seq_along(ph$gtable$grobs)) {
    if (!is.null(ph$gtable$grobs[[i]]$gp$font)) ph$gtable$grobs[[i]]$gp$font <- 1L
    for (j in seq_along(ph$gtable$grobs[[i]]$children)) {
      if (!is.null(ph$gtable$grobs[[i]]$children[[j]]$gp$font)) {
        ph$gtable$grobs[[i]]$children[[j]]$gp$font <- 1L
      }
    }
  }
  ph
}

#* 一张 pheatmap 出 PDF + PNG，尺寸由调用方给
save_pheatmap <- function(ph, stem, width, height) {
  cairo_pdf(paste0(stem, ".pdf"), width = width, height = height, family = plot_font_family)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
  png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 300,
      family = plot_font_family, type = "cairo")
  grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
}

#* 文件名里的时间窗标签：空格与斜杠都不进文件名
slug <- function(x) gsub("[^A-Za-z0-9]+", "-", x)

#* 特征的短标签：多子图排版时全称会把版面吃光，但仍必须是自然语言，不能退回程序字段名
FEATURE_SHORT <- c(
  "Mean Log Effective Population Size" = "Mean Level",
  "Net Change in Log Effective Population Size" = "Net Change",
  "Mean Rate of Change" = "Mean Rate",
  "Fluctuation Amplitude" = "Amplitude",
  "Time of Maximum Expansion" = "Expansion Time",
  "Maximum Expansion Rate" = "Expansion Rate",
  "Time of Maximum Contraction" = "Contraction Time",
  "Maximum Contraction Rate" = "Contraction Rate"
)
short_feature <- function(x) ifelse(x %in% names(FEATURE_SHORT), FEATURE_SHORT[x], x)

#* 打开一个空设备：tidyplots 等包在保存前会隐式画一次，没有当前设备时 R 会在工作目录
#* 里落下一个空白的 Rplots.pdf。开一个 null device 让这些隐式绘图无处可去。
if (is.null(grDevices::dev.list())) grDevices::pdf(NULL)

#* 色条：配置里可以写 diverging / sequential，也可以直接给一串逗号分隔的颜色
resolve_palette <- function(spec, n = 100) {
  spec <- trimws(spec %||% "")
  cols <- if (grepl(",", spec)) {
    trimws(strsplit(spec, ",")[[1]])
  } else if (identical(spec, "sequential")) {
    pal_sequential
  } else if (identical(spec, "diverging") || !nzchar(spec)) {
    pal_diverging
  } else if (exists(spec) && is.character(get(spec))) {
    get(spec)
  } else {
    stop(sprintf("无法解析的色条设置：%s", spec), call. = FALSE)
  }
  colorRampPalette(cols)(n)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

#* 共享图例：一条**真渐变**色条（rasterGrob 插值，不是几十个小矩形拼出来的），
#* 下面接分组色块。所有子图共用一份，既省地方又不会每张子图各画一遍
shared_legend <- function(cols = NULL, lim = NULL, value_title = "Value",
                          groups = NULL, group_cols = NULL,
                          group_title = "Group") {
  gp_txt <- grid::gpar(fontfamily = plot_font_family, fontface = 1, fontsize = 8)
  parts <- list()
  y_top <- 0.96
  if (!is.null(cols) && !is.null(lim)) {
    bar_top <- y_top - 0.05
    bar_h <- 0.34
    parts <- c(parts, list(
      grid::textGrob(value_title, x = 0.06, y = y_top, hjust = 0, gp = gp_txt),
      grid::rasterGrob(rev(cols), x = 0.10, y = bar_top - bar_h / 2, width = 0.16,
                       height = bar_h, hjust = 0, vjust = 0.5, interpolate = TRUE),
      grid::textGrob(sprintf("%.1f", lim), x = 0.30, y = bar_top, hjust = 0, gp = gp_txt),
      grid::textGrob("0", x = 0.30, y = bar_top - bar_h / 2, hjust = 0, gp = gp_txt),
      grid::textGrob(sprintf("%.1f", -lim), x = 0.30, y = bar_top - bar_h,
                     hjust = 0, gp = gp_txt)))
    y_top <- bar_top - bar_h - 0.08
  }
  if (!is.null(groups)) {
    parts <- c(parts, list(grid::textGrob(group_title, x = 0.06, y = y_top,
                                          hjust = 0, gp = gp_txt)))
    for (i in seq_along(groups)) {
      y <- y_top - 0.05 * i
      parts <- c(parts, list(
        grid::rectGrob(x = 0.10, y = y, width = 0.10, height = 0.03, hjust = 0,
                       gp = grid::gpar(fill = group_cols[[groups[i]]], col = NA)),
        grid::textGrob(groups[i], x = 0.24, y = y, hjust = 0, gp = gp_txt)))
    }
  }
  do.call(grid::grobTree, parts)
}

#* 一组子图排成每行 n_col 个的整版图，并按长边限制 PNG 分辨率——子图很多时
#* 300 dpi 会生成上亿像素的图片
save_panels <- function(grobs, stem, n_col = 3, panel_width = 3.6,
                        panel_height = 3.2, title = NULL, legend = NULL,
                        legend_width = 1.6) {
  n_row <- ceiling(length(grobs) / n_col)
  width <- panel_width * min(n_col, length(grobs))
  height <- panel_height * n_row
  combined <- gridExtra::arrangeGrob(grobs = grobs, ncol = n_col, top = title)
  if (!is.null(legend)) {
    # 图例固定高度顶在上方：直接并到整版右侧的话，色条会被整版高度按比例拉成一条巨柱
    legend_height <- min(3.4, height)
    legend <- gridExtra::arrangeGrob(
      legend, grid::nullGrob(), ncol = 1,
      heights = grid::unit(c(legend_height, height - legend_height), "in"))
    combined <- gridExtra::arrangeGrob(
      combined, legend, ncol = 2,
      widths = grid::unit(c(width, legend_width), "in"))
    width <- width + legend_width
  }
  cairo_pdf(paste0(stem, ".pdf"), width = width, height = height,
            family = plot_font_family, onefile = TRUE)
  grid::grid.newpage(); grid::grid.draw(combined); dev.off()
  dpi <- max(60, min(300, 8000 / max(width, height)))
  png(paste0(stem, ".png"), width = width, height = height, units = "in", res = dpi,
      family = plot_font_family, type = "cairo")
  grid::grid.newpage(); grid::grid.draw(combined); dev.off()
}
