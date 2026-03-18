library(ggtree)

tree <- read.tree("/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/output/merge/merged_ml_tree.nwk")

tree  |> ggtree()

#* ==========适用于PAML版本==========
library(tidyverse)
library(ggtree)
library(treeio)

# 设置工作目录
setwd("/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/script")
# 加载树文件
WGS_tree <- read.tree("/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/output/merge/merged_ml_tree.nwk")
# 文件路径可按需修改
group_file      <- "./conf/group.txt"                                #! 第一列是label(ID)，第二列是group(单倍群)，不包含表头
color_file      <- "./conf/color.txt"                                #! 第一列是group(单倍群)，第二列是color，不包含表头

##* =========================分支着色(如果没有分组请删除下列的aes(color = group))=========================
##* =========================分支着色(如果没有分组请删除下列的scale_color_manual(values = group_colors) =========================
# 读取分组信息
group_df        <- read.table(group_file, 
                              header = FALSE, 
                              sep = "\t", 
                              stringsAsFactors = FALSE, 
                              col.names = c("label", "haplo"))
group_list      <- split(group_df$label, group_df$haplo)              # 将分组信息转换为列表
# 读取颜色信息
color_df        <- read.table(color_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = '', col.names = c("group", "color"))
group_colors    <- setNames(color_df$color, color_df$group)           # 将颜色信息转换为命名向量

# 读取树并分组
tree_grouped <- groupOTU(WGS_tree, group_list)

p0 <- ggtree(
    tree_grouped,                                                     #* 系统发育树对象
    aes(color = group),                                               #* 分支颜色
    layout = 'rectangular',                                           #* 布局方式
    open.angle = 5,                                                   #* 开角度
    lwd = 0.2                                                         #* branch宽度
) +
scale_color_manual(values = group_colors) +                           #* 如果没有分支着色请删除此行
geom_tiplab(
    size = 0.5,                                                       #* tiplab文字大小
    align = TRUE,                                                     #* 是否对齐
    hjust = -0.5                                                      #* 水平调整
) 
p0
##* 保存
ggsave(p0, filename = '系统发育树.pdf', units = "cm", height = 20, width = 15 )


#* ==========适用于BEAST MCCtree版本==========

library(ggplot2)                                                      #* 绘图基础库
library(ggtree)                                                       #* 系统发育树可视化
library(treeio)                                                       #* BEAST树数据读取

base_dir <- "/mnt/c/Users/Administrator/Desktop"                      #* 工作目录基准
tree_file <- file.path(base_dir, "基线-CloneML_rename.tree")            #* 树文件路径
group_file <- file.path(base_dir, "conf/group.txt")                   #* 分组文件路径
color_file <- file.path(base_dir, "conf/color.txt")                   #* 颜色文件路径

read_tsv <- function(path, cols) {
    read.table(
        path,
        header = FALSE,
        sep = "\t",
        stringsAsFactors = FALSE,
        col.names = cols,
        comment.char = ""
    )                                                                 #* 读取TSV
}

group_df <- read_tsv(group_file, c("label", "group"))                 #* 加载分组信息
group_list <- split(group_df$label, group_df$group)                   #* 分组转为列表
color_df <- read_tsv(color_file, c("group", "color"))                 #* 加载颜色信息
group_colors <- setNames(color_df$color, color_df$group)              #* 生成颜色映射

beast_tree <- read.beast(tree_file)                                   #* 读取BEAST树
tree_grouped <- groupOTU(beast_tree, group_list)                      #* 按分组着色分支
tree_df <- ggtree::fortify(tree_grouped)                              #* 提取节点坐标与注释
max_height <- max(tree_df$x, na.rm = TRUE)                            #* 树高度最大值
tree_df$node_to_tip_distance <- max_height - tree_df$x                #* 节点到现生距离

hpd_bounds <- t(vapply(tree_df$height_0.95_HPD, function(x) {
    if (length(x) == 2 && all(is.finite(x))) x else c(NA_real_, NA_real_)
}, numeric(2)))
tree_df$hpd_lower <- hpd_bounds[, 1]                                  #* HPD下限
tree_df$hpd_upper <- hpd_bounds[, 2]                                  #* HPD上限
tree_df$hpd_xmin <- max_height - tree_df$hpd_upper                    #* HPD映射到绘图坐标下限
tree_df$hpd_xmax <- max_height - tree_df$hpd_lower                    #* HPD映射到绘图坐标上限
tree_df$hpd_label <- ifelse(
    is.finite(tree_df$hpd_lower) & is.finite(tree_df$hpd_upper),
    sprintf("[%.0f, %.0f]", tree_df$hpd_lower, tree_df$hpd_upper),
    NA_character_
)                                                                     #* HPD上下限标签文本

break_5k <- 5000                                                      #* 5k步长
break_10k <- 10000                                                    #* 10k步长
time_breaks <- function()                                             #* x轴刻度
    max_height - seq(0, max_height, by = break_5k)
time_labels <- function(x)                                            #* x轴标签
    format(round(max_height - x, 0), big.mark = ",")
vline_5k <- (                                                         #* 5000年网格
    max_height - seq(0, max_height, by = break_5k)
)
vline_10k <- (                                                        #* 10000年网格
    max_height - seq(0, max_height, by = break_10k)
)

p0 <- ggtree(                                                         #* 初始化树图
    tree_grouped,                                                     #* 树对象
    aes(color = group),                                               #* 分支按分组着色
    layout = "rectangular",                                           #* 矩形布局
    open.angle = 5,                                                   #* 开角度
    lwd = 0.2                                                         #* 分支线宽
) +                                                                   #* 继续叠加图层
scale_color_manual(values = group_colors) +                           #* 指定分组颜色
geom_errorbar(                                                        #* 绘制HPD条
    data = subset(tree_df, !isTip),                                   #* 仅内部节点
    aes(y = y, xmin = hpd_xmin, xmax = hpd_xmax),                     #* HPD范围
    orientation = "y",                                                #* 水平误差条
    width = 0.2,                                                      #* 误差条高度
    color = "#9fd3f0",                                                #* 误差条颜色
    alpha = 0.8,                                                      #* 透明度
    linewidth = 0.25                                                  #* 线宽
) +                                                                   #* 继续叠加图层
geom_tiplab(                                                          #* 末端标签
    size = 0.5,                                                       #* 字号
    align = TRUE,                                                     #* 对齐标签
    hjust = 0,                                                        #* 水平对齐
    offset = max_height * 0.02                                        #* 右移偏移
) +                                                                   #* 继续叠加图层
geom_nodepoint(size = 0.3, alpha = 0.5) +                             #* 节点点位
geom_nodelab(                                                         #* 分支长度标签
    data = subset(tree_df, !isTip),                                   #* 仅内部节点
    aes(label = round(branch.length, 3)),                             #* 分支长度文本
    hjust = -0.3,                                                     #* 水平偏移
    size = 0.7                                                        #* 字号
) +                                                                   #* 继续叠加图层
geom_nodelab(                                                         #* 到现生距离标签
    data = subset(tree_df, !isTip),                                   #* 仅内部节点
    aes(label = round(node_to_tip_distance, 3)),                      #* 距离文本
    hjust = 1,                                                        #* 水平偏移
    color = "#000000",                                                #* 文字颜色
    size = 0.7                                                        #* 字号
) +                                                                   #* 继续叠加图层
geom_nodelab(                                                         #* HPD上下限标签
    data = subset(tree_df, !isTip),                                   #* 仅内部节点
    aes(label = hpd_label),                                           #* HPD文本
    hjust = 1.2,                                                      #* 水平偏移
    vjust = -0.4,                                                     #* 垂直偏移
    color = "#3b3b3b",                                                #* 文字颜色
    size = 0.6                                                        #* 字号
) +                                                                   #* 继续叠加图层
theme(                                                                #* 主题设置
    axis.text.x = element_text(size = 2),                             #* x轴刻度字号
    axis.title.x = element_text(size = 2),                            #* x轴标题字号
    axis.text.y = element_text(size = 2),                             #* y轴刻度字号
    axis.title.y = element_text(size = 2),                            #* y轴标题字号
    plot.margin = margin(5, 20, 5, 5),                                #* 绘图区边距
    legend.position = "bottom"                                        #* 图例位置
) +                                                                   #* 继续叠加图层
scale_x_continuous(                                                   #* x轴刻度
    expand = expansion(mult = c(0.02, 0.12)),                         #* x轴留白
    breaks = function(limits) time_breaks(),                          #* x轴刻度位置
    labels = time_labels,                                             #* x轴刻度文本
    name = "Time (years ago)"                                         #* x轴标题
) +                                                                   #* 继续叠加图层
scale_y_continuous(                                                   #* y轴刻度
    expand = expansion(mult = c(0.01, 0.01))                          #* y轴留白
) +                                                                   #* 继续叠加图层
coord_cartesian(clip = "off") +                                       #* 不裁剪标签
geom_vline(                                                           #* 5000年网格线
    xintercept = vline_5k,                                            #* 位置
    color = "#1c9cc6",                                                #* 颜色
    linetype = "solid",                                               #* 线型
    alpha = 0.3,                                                      #* 透明度
    linewidth = 0.05                                                  #* 线宽
) +                                                                   #* 继续叠加图层
geom_vline(                                                           #* 10000年网格线
    xintercept = vline_10k,                                           #* 位置
    color = "#196c88",                                                #* 颜色
    linetype = "solid",                                               #* 线型
    alpha = 0.5,                                                      #* 透明度
    linewidth = 0.1                                                   #* 线宽
)                                                                     #* 图层结束

ggsave(p0, filename = file.path(base_dir, "系统发育树.pdf"), units = "cm", height = 20, width = 15) #* 导出PDF

library(ggtreeExtra)
#* 把单倍群颜色环带添加到第一层热图
p1 <- p0 +
    new_scale_fill() +
    geom_fruit(
        data = group_df,                            #* 使用分组数据框
        geom = geom_tile,                           #* geom_tile 用于绘制热图
        mapping = aes(y = label, fill = haplo),     #* 映射标签和单倍群
        offset = 0.001,                             #* 热图偏移量
        width = 0.0005                              #* 热图宽度
    ) +
    geom_fruit(
        data = group_df |> distinct(haplo, .keep_all = TRUE), #* 获取唯一的单倍群标签
        geom = geom_text,                                     #* 使用文本几何对象 
        mapping = aes(y = label, label = haplo),              #* 映射标签和单倍群
        size = 2,                                             #* 字体大小
        color = "#ffffff",                                  #* 字体颜色
        offset = 0.001,                                       #* 字体偏移量
        pwidth = 0.0005                                       #* 字体宽度
    ) +
    scale_fill_manual(
        values = group_colors                                 #* 使用预定义的颜色映射
    )
p1
ggsave(p1, filename = '../output/系统发育树.pdf', units = "cm", height = 20, width = 15 )
