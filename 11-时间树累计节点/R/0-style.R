#* =====字体注册=====
# 直接注册 R 自带的 Arial 字体度量，使 grDevices::pdf() 输出文字可编辑的 Arial。
# 系统没有这些度量文件时回落到默认字体，不让整条流程因字体而失败。
afm_dir <- file.path(R.home("library"), "grDevices", "afm")
afm_files <- file.path(afm_dir, c("ArialMT.afm.gz", "ArialMT-Bold.afm.gz",
                                  "ArialMT-Italic.afm.gz", "ArialMT-BoldItalic.afm.gz"))
plot_family <- if (all(file.exists(afm_files))) {
  grDevices::pdfFonts(Arial = grDevices::Type1Font("Arial", afm_files))
  "Arial"
} else {
  message("未找到 Arial 字体度量，回落到默认字体")
  ""
}

#* =====主题=====
# 复用 tidyplots 的 theme_tidyplot 口径：白底、只留 x/y 轴线、全部文字同一字号同一颜色。
# tidyplots 的主题函数只接受它自己构造的图对象，故这里按同样的元素值直接组装。
theme_tidyplot_like <- function(fontsize = 7, family = "") {
  text_element <- element_text(size = fontsize, family = family, face = "plain", colour = "black")
  theme_grey(base_size = fontsize, base_family = family) +
    theme(
      plot.margin = unit(c(1, 1, 1, 1), "mm"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      strip.background = element_rect(fill = NA, colour = NA),
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      text = text_element,
      axis.text = text_element,
      axis.title = text_element,
      strip.text = text_element,
      legend.title = text_element,
      legend.text = text_element,
      plot.title = element_text(size = fontsize, family = family, face = "plain",
                                colour = "black", hjust = 0.5, vjust = 0.5),
      legend.position = "none"
    )
}

#* =====命令行参数=====
# 用法：get_arg("--ridge-scale")；numeric = TRUE 时转成数值，split 非空时按该分隔符拆成向量
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, numeric = FALSE, split = NULL) {
  value <- args[which(args == name) + 1]
  if (!is.null(split)) value <- str_split_1(value, split)
  if (numeric) as.numeric(value) else value
}

#* =====展示名=====
# 文件名用连字符，不落任何下划线命名
slugify <- function(x) str_replace_all(str_replace_all(x, "[^A-Za-z0-9]+", "-"), "^-|-$", "")
