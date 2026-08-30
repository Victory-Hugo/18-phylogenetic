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
#* lim 给一个数是对称色标 [-lim, lim]，给两个数就是 [lo, hi]——同簇频率这类
#* 取值在 0-1 的量没有"以 0 为中心"的读法，硬套对称色条会读错
shared_legend <- function(cols = NULL, lim = NULL, value_title = "Value",
                          groups = NULL, group_cols = NULL,
                          group_title = "Group") {
  gp_txt <- grid::gpar(fontfamily = plot_font_family, fontface = 1, fontsize = 8)
  parts <- list()
  y_top <- 0.96
  if (!is.null(cols) && !is.null(lim)) {
    rng <- if (length(lim) >= 2) sort(as.numeric(lim[1:2])) else c(-abs(lim), abs(lim))
    digits <- if (diff(rng) <= 3) 2 else 1
    lab <- function(v) formatC(v, format = "f", digits = digits)
    bar_top <- y_top - 0.05
    bar_h <- 0.34
    parts <- c(parts, list(
      grid::textGrob(value_title, x = 0.06, y = y_top, hjust = 0, gp = gp_txt),
      grid::rasterGrob(rev(cols), x = 0.10, y = bar_top - bar_h / 2, width = 0.16,
                       height = bar_h, hjust = 0, vjust = 0.5, interpolate = TRUE),
      grid::textGrob(lab(rng[2]), x = 0.30, y = bar_top, hjust = 0, gp = gp_txt),
      grid::textGrob(lab(mean(rng)), x = 0.30, y = bar_top - bar_h / 2, hjust = 0,
                     gp = gp_txt),
      grid::textGrob(lab(rng[1]), x = 0.30, y = bar_top - bar_h,
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

#* =====A4 纵向分页=====
A4_WIDTH <- 8.27
A4_HEIGHT <- 11.69

#* 内容比 A4 宽时不减列、不换方向，而是整页等比放大出图、再统一缩回 A4：
#* 减列会把 40 个时间段摊成几十页，直接截断则打印时丢内容
shrink_pdf_to_a4 <- function(path) {
  if (nzchar(Sys.which("pdfjam"))) {
    tmp <- paste0(path, ".a4.pdf")
    ok <- system2("pdfjam", c("--paper", "a4paper", "--quiet", "--outfile",
                              shQuote(tmp), shQuote(path)),
                  stdout = FALSE, stderr = FALSE)
    if (ok == 0 && file.exists(tmp)) {
      file.rename(tmp, path)
      return(TRUE)
    }
    unlink(tmp)
  }
  if (nzchar(Sys.which("gs"))) {
    tmp <- paste0(path, ".a4.pdf")
    ok <- system2("gs", c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
                          "-sPAPERSIZE=a4", "-dFIXEDMEDIA", "-dPDFFitPage",
                          paste0("-sOutputFile=", tmp), path),
                  stdout = FALSE, stderr = FALSE)
    if (ok == 0 && file.exists(tmp)) {
      file.rename(tmp, path)
      return(TRUE)
    }
    unlink(tmp)
  }
  FALSE
}

#* 一组子图排成每行 n_col 个的整版图，按 A4 纵向分页：PDF 是多页的，PNG 每页一个文件
#* （只有一页时文件名不带页码，与旧输出保持一致）。长边限制 PNG 分辨率——子图很多时
#* 300 dpi 会生成上亿像素的图片
save_panels <- function(grobs, stem, n_col = 3, panel_width = 3.6,
                        panel_height = 3.2, title = NULL, legend = NULL,
                        legend_width = 1.6) {
  if (!length(grobs)) stop("save_panels 收到空的子图列表", call. = FALSE)
  n_used_col <- min(n_col, length(grobs))
  body_width <- panel_width * n_used_col
  page_width <- body_width + if (is.null(legend)) 0 else legend_width
  #* 页面先按 A4 的宽高比放大到内容宽度，最后整份 PDF 再等比缩回真正的 A4
  page_scale <- max(1, page_width / A4_WIDTH)
  page_width <- A4_WIDTH * page_scale
  page_height <- A4_HEIGHT * page_scale
  title_h <- if (is.null(title)) 0 else 0.5 * page_scale
  foot_h <- 0.3 * page_scale
  rows_per_page <- max(1, floor((page_height - title_h - foot_h) / panel_height))
  per_page <- rows_per_page * n_used_col
  n_page <- ceiling(length(grobs) / per_page)
  body_h <- rows_per_page * panel_height

  page_grob <- function(k) {
    idx <- seq((k - 1) * per_page + 1, min(k * per_page, length(grobs)))
    g <- grobs[idx]
    if (length(g) < per_page) {
      g <- c(g, replicate(per_page - length(g), grid::nullGrob(), simplify = FALSE))
    }
    body <- gridExtra::arrangeGrob(
      grobs = g, ncol = n_used_col,
      widths = grid::unit(rep(panel_width, n_used_col), "in"),
      heights = grid::unit(rep(panel_height, rows_per_page), "in"))
    if (!is.null(legend)) {
      # 图例固定高度顶在上方：直接并到整版右侧的话，色条会被整版高度按比例拉成一条巨柱
      legend_height <- min(3.4 * page_scale, body_h)
      leg <- gridExtra::arrangeGrob(
        legend, grid::nullGrob(), ncol = 1,
        heights = grid::unit(c(legend_height, body_h - legend_height), "in"))
      body <- gridExtra::arrangeGrob(
        body, leg, ncol = 2,
        widths = grid::unit(c(body_width, legend_width), "in"))
    }
    foot <- if (n_page > 1) {
      grid::textGrob(sprintf("Page %d of %d", k, n_page), y = 0.5,
                     gp = grid::gpar(fontfamily = plot_font_family, fontface = 1,
                                     fontsize = 8, col = "grey40"))
    } else {
      grid::nullGrob()
    }
    parts <- list(if (is.null(title)) grid::nullGrob() else title, body, foot,
                  grid::nullGrob())
    gridExtra::arrangeGrob(
      grobs = parts, ncol = 1,
      heights = grid::unit(c(title_h, body_h, foot_h,
                             max(page_height - title_h - body_h - foot_h, 0)), "in"))
  }

  #* 上一次运行留下的分页 PNG 必须先清掉：页数变少时，旧页会冒充本次的输出
  unlink(c(paste0(stem, ".png"), Sys.glob(sprintf("%s-p*.png", stem))))

  pdf_path <- paste0(stem, ".pdf")
  cairo_pdf(pdf_path, width = page_width, height = page_height,
            family = plot_font_family, onefile = TRUE)
  for (k in seq_len(n_page)) {
    grid::grid.newpage(); grid::grid.draw(page_grob(k))
  }
  dev.off()
  if (page_scale > 1 && !shrink_pdf_to_a4(pdf_path)) {
    message(sprintf("找不到 pdfjam / gs，%s 的页面保持 %.2f x %.2f in（A4 的 %.2f 倍）",
                    basename(pdf_path), page_width, page_height, page_scale))
  }

  dpi <- max(60, min(300, 8000 / max(page_width, page_height)))
  for (k in seq_len(n_page)) {
    png_path <- if (n_page == 1) paste0(stem, ".png") else
      sprintf("%s-p%02d.png", stem, k)
    png(png_path, width = page_width, height = page_height, units = "in", res = dpi,
        family = plot_font_family, type = "cairo")
    grid::grid.newpage(); grid::grid.draw(page_grob(k))
    dev.off()
  }
  invisible(n_page)
}
