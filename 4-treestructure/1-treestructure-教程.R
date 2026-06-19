# ================================
# treestructure 示例分析脚本
# 功能：读取示例时间树，识别隐藏群体结构，并输出结果表
# ================================

# 安装 treestructure
# 如果已经安装过，可以不用重复运行
# install.packages("devtools")
# devtools::install_github("emvolz-phylodynamics/treestructure")

# 加载系统发育树处理包
library(ape)

# 加载 treestructure 包
library(treestructure)

# 设置工作目录
# 后续输出文件会保存在这个目录下
setwd('/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/4-treestructure/')

# 读取 treestructure 包自带的模拟树
# sim.nwk 是一个 Newick 格式的示例树文件
tree <- read.tree(system.file("sim.nwk", package = "treestructure"))

# 运行 treestructure
# trestruct() 会根据时间树中的共祖模式识别隐藏群体结构
s <- trestruct(tree)

# 打印分析摘要
# 包括 cluster 数量、partition 数量，以及每类包含的样本数
print(s)

# 将 treestructure 结果转换为数据框
# result 中包含每个样本对应的 cluster 和 partition
result <- as.data.frame(s)

# 查看每个 cluster 中的样本数量
table(result$cluster)

# 查看每个 partition 中的样本数量
table(result$partition)

# 输出 TSV 文件
# sep = "\t" 表示使用制表符分隔
# quote = FALSE 表示不为字符加引号
# row.names = FALSE 表示不输出行名
write.table(
  result,
  file = "treestructure_result.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 绘制 treestructure 分群结果
plot(s)
