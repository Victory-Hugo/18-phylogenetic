library(tidyverse)
library(ggtree)
library(treeio)

tree_path <- "/mnt/l/6_起源地混合地/8-第一轮修稿/5-拆树PAML合树/data/VFT_rooted.tree"
tree <- read.tree(tree_path)
tree  |> ggtree()

# 设置工作目录
setwd("/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/script")
# 加载树文件
WGS_tree <- read.tree(tree_path)                                              #* 读取树文件，路径可按需修改
# 文件路径可按需修改
group_file      <- "./conf/haplogrep_group.txt"                                #! 第一列是label(ID)，第二列是group(单倍群)，不包含表头
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
