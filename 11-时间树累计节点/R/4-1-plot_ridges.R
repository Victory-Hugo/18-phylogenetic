library(tidyverse)
library(ggridges)
library(scales)
library(patchwork)

source("R/0-style.R")

#* =====参数=====
ridge_density_path  <- get_arg("--ridge-density")
group_design_path   <- get_arg("--group-design")
output_figure       <- get_arg("--output-figure")
fill_mode           <- get_arg("--fill-mode")
colour_ramp         <- get_arg("--colour-ramp", split = ",")
ridge_scale         <- get_arg("--ridge-scale", numeric = TRUE)
ridge_alpha         <- get_arg("--ridge-alpha", numeric = TRUE)
outline_width       <- get_arg("--outline-width", numeric = TRUE)
min_height_fraction <- get_arg("--min-height-fraction", numeric = TRUE)
baseline_colour     <- get_arg("--baseline-colour")
outline_colour      <- get_arg("--outline-colour")
reference_colour    <- get_arg("--reference-colour")
interval_colour     <- get_arg("--interval-colour")
interval_alpha      <- get_arg("--interval-alpha", numeric = TRUE)
reference_recent    <- get_arg("--reference-lines-recent", numeric = TRUE, split = ",")
reference_full      <- get_arg("--reference-lines-full", numeric = TRUE, split = ",")
font_size           <- get_arg("--font-size", numeric = TRUE)
panel_mm            <- get_arg("--panel-mm", numeric = TRUE)
label_mm            <- get_arg("--label-mm", numeric = TRUE)
min_row_mm          <- get_arg("--min-row-mm", numeric = TRUE)

theme_set(theme_tidyplot_like(fontsize = font_size, family = plot_family))
dir.create(output_figure, recursive = TRUE, showWarnings = FALSE)

#* =====读取密度长表=====
df1 <- read_tsv(ridge_density_path, show_col_types = FALSE)
df2 <- read_tsv(group_design_path, show_col_types = FALSE) %>%
  select(Grouping, Group, Order, Colour, Outline)
df3 <- inner_join(df1, df2, by = c("Grouping", "Group")) %>%
  mutate(x_max = as.numeric(str_extract(`Time window`, "[0-9.]+")))

groupings <- unique(df2$Grouping)
windows <- df3 %>% distinct(`Time window`, x_max) %>% arrange(x_max)

#* =====单个面板=====
# 面板强制正方形（aspect.ratio = 1）；同一维度的各类别共用一个高度换算系数，脊线高度真正可比。
build_panel <- function(df, x_max) {
  levels_tbl <- df %>%
    distinct(Group, Order, Samples, `Share of events below the resolution limit`) %>%
    arrange(desc(Order)) %>%
    mutate(y_pos = row_number(),
           axis_label = paste0(str_wrap(Group, width = 22), "\nn = ", comma(Samples),
                               if_else(`Share of events below the resolution limit` > 0,
                                       paste0(", ", percent(
                                         `Share of events below the resolution limit`,
                                         accuracy = 1), " excluded"), "")))
  n_rows <- nrow(levels_tbl)

  df4 <- df %>%
    left_join(select(levels_tbl, Group, y_pos), by = "Group") %>%
    arrange(desc(y_pos), `Node age (kya)`)
  effective_scale <- min(ridge_scale, 0.9 * n_rows)
  height_factor <- effective_scale /
    max(df4$`Ridge height upper`, df4$`Ridge height`, na.rm = TRUE)
  # 低于全图最高脊线该比例的部分不画：在共用高度尺子下这部分本来就细到看不见，
  # 裁掉之后脊线左端收成尖角，而不是一路平铺成长条
  ridge_min_height <- min_height_fraction * effective_scale

  p <- ggplot() +
    geom_hline(yintercept = seq_len(n_rows), linewidth = 0.25, colour = baseline_colour) +
    geom_vline(xintercept = if (x_max <= 10) reference_recent else reference_full,
               colour = reference_colour, linewidth = 0.2, linetype = "dotted")

  # 抽稀重复之间的区间垫在脊线之下，取该条脊线自己的填充色再调淡，不用中性灰
  if (!all(is.na(df$`Ridge height lower`))) {
    band <- filter(df4, `Ridge height upper` * height_factor >= ridge_min_height)
    band_aes <- aes(x = `Node age (kya)`,
                    ymin = y_pos + `Ridge height lower` * height_factor,
                    ymax = y_pos + `Ridge height upper` * height_factor, group = Group)
    p <- p + if (fill_mode == "category") {
      # age 模式的填充色由年龄决定，没有属于某条脊线的单色，那里才回落到中性色
      geom_ribbon(data = band, modifyList(band_aes, aes(fill = Colour)),
                  alpha = interval_alpha)
    } else {
      geom_ribbon(data = band, band_aes, fill = interval_colour, alpha = interval_alpha)
    }
  }

  # 填充与描边分开画：填充带透明度，描边单独一条不透明的线，因此描边始终比填充更深更实
  outline_data <- filter(df4, `Ridge height` * height_factor >= ridge_min_height)
  if (fill_mode == "category") {
    # 每条脊线一个纯色，颜色沿色带按纵轴顺序插值；描边取同色压暗一档
    p <- p +
      geom_ridgeline(
        data = df4,
        aes(x = `Node age (kya)`, y = y_pos, height = `Ridge height` * height_factor,
            group = Group, fill = Colour),
        colour = NA, alpha = ridge_alpha, min_height = ridge_min_height) +
      geom_line(data = outline_data,
                aes(x = `Node age (kya)`, y = y_pos + `Ridge height` * height_factor,
                    group = Group, colour = Outline), linewidth = outline_width) +
      scale_fill_identity() +
      scale_colour_identity()
  } else {
    # 脊线内部沿节点年龄连续渐变，色带范围即当前时间窗口
    p <- p +
      geom_ridgeline_gradient(
        data = df4,
        aes(x = `Node age (kya)`, y = y_pos, height = `Ridge height` * height_factor,
            group = Group, fill = `Node age (kya)`),
        colour = NA, alpha = ridge_alpha, min_height = ridge_min_height) +
      geom_line(data = outline_data,
                aes(x = `Node age (kya)`, y = y_pos + `Ridge height` * height_factor,
                    group = Group), colour = outline_colour, linewidth = outline_width) +
      scale_fill_gradientn(colours = colour_ramp, limits = c(0, x_max),
                           oob = scales::squish, guide = "none")
  }

  # 高度标尺：面板左上角的竖直参照线段。同一面板内各类别共用一个高度换算系数，
  # 脊线高度本身可比，这把标尺给出它的物理单位
  bar_value <- signif(max(df4$`Ridge height`, na.rm = TRUE) / 4, 1)
  bar_height <- bar_value * height_factor
  top_row <- filter(df4, y_pos == n_rows)
  top_row_height <- max(c(top_row$`Ridge height`, top_row$`Ridge height upper`),
                        na.rm = TRUE) * height_factor
  left_height <- top_row %>%
    filter(`Node age (kya)` >= x_max * 0.55) %>%
    summarise(v = max(c(`Ridge height`, `Ridge height upper`), na.rm = TRUE)) %>%
    pull(v) * height_factor
  bar_x <- x_max * 0.97
  bar_y <- n_rows + left_height + 0.06
  # 纵轴范围已被最上一条脊线自身的高度撑开，这里只需补出标尺高过该脊线的那一点点
  top_pad <- max(0, bar_y + bar_height + 0.16 - (n_rows + top_row_height)) + 0.1

  p <- p +
    annotate("segment", x = bar_x, xend = bar_x, y = bar_y, yend = bar_y + bar_height,
             colour = "black", linewidth = 0.35) +
    annotate("text", x = bar_x - x_max * 0.025, y = bar_y + bar_height / 2,
             label = comma(bar_value), family = plot_family, fontface = "plain",
             colour = "black", size = font_size / .pt, hjust = 0, vjust = 0.5) +
    scale_x_reverse(breaks = pretty(c(0, x_max), n = 6),
                    expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(breaks = levels_tbl$y_pos,
                       labels = levels_tbl$axis_label,
                       expand = expansion(add = c(0.1, top_pad))) +
    coord_cartesian(xlim = c(x_max, 0), clip = "off") +
    labs(x = "Node age (kya)", y = NULL) +
    theme(aspect.ratio = 1,
          axis.text.y = element_text(hjust = 1, vjust = 0.5, lineheight = 1.05),
          axis.ticks.y = element_blank())
}

#* =====每个维度一张图：各时间窗口横排为子图=====
# standalone = FALSE 时把维度名写进每个面板的标题：patchwork 嵌套会丢掉子图自己的总标题，
# 合成图只能靠面板标题标出维度。
build_figure <- function(grouping, standalone = TRUE) {
  panels <- list()
  for (i in seq_len(nrow(windows))) {
    df <- filter(df3, Grouping == grouping, `Time window` == windows$`Time window`[i])
    if (nrow(df) == 0) next
    panels[[length(panels) + 1]] <- build_panel(df, windows$x_max[i]) +
      labs(title = if (standalone) windows$`Time window`[i]
           else paste0(grouping, ", ", str_to_lower(windows$`Time window`[i])))
  }
  if (length(panels) == 0) return(NULL)
  unit_label <- str_wrap(first(filter(df3, Grouping == grouping)$`Ridge height unit`), width = 30)
  panels <- map(panels, ~ .x + labs(y = unit_label) +
                  theme(axis.title.y = element_text(angle = 90, vjust = 1)))
  p <- wrap_plots(panels, nrow = 1)
  if (standalone) {
    p <- p + plot_annotation(title = paste0("Coalescent event density by ", grouping),
                             theme = theme(plot.title = element_text(hjust = 0, lineheight = 1.2)))
  }
  p
}

save_figure <- function(plot, stem, width_mm, height_mm, dpi = 300) {
  grDevices::pdf(file.path(output_figure, paste0(stem, ".pdf")),
                 width = width_mm / 25.4, height = height_mm / 25.4,
                 family = plot_family, useDingbats = FALSE)
  print(plot)
  dev.off()
  ggsave(file.path(output_figure, paste0(stem, ".png")), plot, width = width_mm,
         height = height_mm, units = "mm", dpi = dpi, limitsize = FALSE)
  message(sprintf("%-52s (%.0f x %.0f mm)", stem, width_mm, height_mm))
}

# 每个维度的类别数不同；面板是正方形，边长按类别数放大，否则类别多时每行太薄、标签会叠在一起
n_groups <- df2 %>% count(Grouping) %>% deframe()
panel_size <- function(rows) max(panel_mm, rows * min_row_mm)

for (i in seq_along(groupings)) {
  p <- build_figure(groupings[i])
  if (is.null(p)) next
  side <- panel_size(n_groups[[groupings[i]]])
  save_figure(p, paste0("⭐4-", i, "-", slugify(groupings[i])),
              nrow(windows) * (side + label_mm), side + 20)
}

#* =====合成图：全部维度纵向排列，与上面同样的横排子图=====
# 各维度的脊线高度换算系数仍各自独立，因此跨维度不宜比高度，每个面板都保留自己的高度标尺。
rows <- keep(map(groupings, build_figure, standalone = FALSE), ~ !is.null(.x))
if (length(rows) > 1) {
  # 各块共用一个边长，取类别最多的那块的需求，否则它会被压扁
  side <- panel_size(max(n_groups))
  save_figure(wrap_plots(rows, ncol = 1), "⭐4-composite-All-groupings",
              nrow(windows) * (side + label_mm), length(rows) * (side + 14) + 8, dpi = 200)
}
