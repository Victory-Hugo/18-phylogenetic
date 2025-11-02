# pplacer 可视化脚本使用指南

## 功能说明

这个脚本用于对 pplacer 的 jplace 文件进行可视化处理，生成：
- 系统发育树的圆形树图（带放置数据）
- 似然权重比 (Likelihood Weight Ratio, LWR) 分布直方图
- 序列放置不确定性展示
- 分支权重可视化
- 增强树文件（Newick 格式）

## 快速开始

### 基础用法

```bash
# 方式 1：直接调用 R 脚本（jplace 文件所在目录为输出目录）
Rscript 2-可视化.R /path/to/file.jplace

# 方式 2：指定输出目录
Rscript 2-可视化.R /mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/output/placer.aln.snp-sites.raxml.jplace "" /mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/output/fig/

# 方式 3：使用 shell 包装脚本
bash run_visualization.sh /path/to/file.jplace
```

### 完整用法

```bash
# 不指定分类文件，输出到输入文件所在目录
Rscript 2-可视化.R data/placer.jplace

# 指定输出目录
Rscript 2-可视化.R data/placer.jplace "" results/

# 指定分类文件和输出目录
Rscript 2-可视化.R data/placer.jplace data/classification.csv results/

# 使用 shell 脚本（推荐）
bash run_visualization.sh data/placer.jplace data/classification.csv results/
```

## 参数说明

### 位置参数

| 参数 | 位置 | 必需 | 说明 |
|-----|------|------|------|
| jplace文件 | 1 | ✓ | pplacer 生成的 jplace 文件路径 |
| 分类文件 | 2 | ✗ | 分类信息 CSV 文件（可选，格式见下文） |
| 输出目录 | 3 | ✗ | 输出文件夹路径，不存在会自动创建（默认为输入文件所在目录） |

### 分类文件格式

分类 CSV 文件应包含两列：
- 第一列：序列标签（需与 jplace 中的 tip labels 对应）
- 第二列：分类信息（group/taxon_group 等）

示例（`classification.csv`）：
```
label,group
seq_001,Proteobacteria
seq_002,Firmicutes
seq_003,Proteobacteria
```

列名可以是以下任意一个（自动识别）：
- label、tip、tip_label、tip.label、taxon、name
- group、Group、taxon_group、category、clade

## 输出文件

脚本会在输出目录生成以下文件：

### 必生成文件
- `jplace_overview.csv` - jplace 文件元信息
- `jplace_fields.csv` - 字段定义
- `df_pplacer_likelihood.csv` - 所有放置的似然信息
- `likelihood_weight_ratio_hist.png` - LWR 分布直方图

### 条件生成文件

#### 圆形树（两个版本）
- `placement_circular_tree_ignore_branch_length.png` - 不考虑分支长度
- `placement_circular_tree_with_branch_length.png` - 考虑分支长度
- `placement_circular_tree_groups_*.png` - 带分类注释的版本（若提供分类文件）

#### 序列放置可视化
- `placement_sequence_<序列名>_*.png` - 各序列的放置不确定性展示

#### 分支权重可视化
- `placement_branch_weights_*.png` - 分支权重图

#### 增强树
- `placement_augmented_tree.nwk` - 包含序列放置的增强树
- `placement_augmented_mapping.csv` - 原始序列名到最终标签的映射

#### 其他可选输出
- `jplace_metadata.csv` - 元数据（若存在）
- `jplace_placement_names.csv` - 序列名信息
- `jplace_placement_candidates.csv` - 所有候选放置
- `jplace_placement_extras.csv` - 额外字段信息

## 使用示例

### 示例 1：基础分析

```bash
# 最简单的用法
Rscript 2-可视化.R my_data/placer.jplace
```

### 示例 2：指定输出目录

```bash
# 输出到 results 文件夹
Rscript 2-可视化.R my_data/placer.jplace "" results/
```

### 示例 3：完整分析（包含分类信息）

```bash
# 提供分类信息，生成带分类注释的树图
Rscript 2-可视化.R my_data/placer.jplace metadata/taxonomy.csv results/
```

### 示例 4：使用 shell 脚本

```bash
# 最推荐的方式，带详细日志
bash run_visualization.sh my_data/placer.jplace metadata/taxonomy.csv results/
```

## 所需 R 包

脚本会自动加载以下包（需提前安装）：

```r
install.packages(c("jsonlite", "ape", "dplyr", "ggplot2", "colorspace", "tibble"))
devtools::install_github("YuLab-SMU/treeio")
devtools::install_github("YuLab-SMU/ggtree")
devtools::install_github("YuLab-SMU/ggtreeExtra")
```

## 故障排除

### 问题 1：`Rscript: command not found`
**解决**：安装 R 或检查 PATH 环境变量

### 问题 2：缺少 R 包
**解决**：在 R 中执行安装命令，见上方"所需 R 包"

### 问题 3：输入文件路径错误
**解决**：确保 jplace 文件路径正确，使用绝对路径或相对路径均可

### 问题 4：输出目录权限不足
**解决**：确保输出目录有写权限，或指定其他位置

## 性能提示

- 较大的 jplace 文件（> 100 MB）处理时间可能较长
- 序列数过多可能导致圆形树难以阅读
- 建议使用 SSD 存储以加快 I/O

## 技术细节

### 默认配置
- 色彩方案：Sunset (colorspace 包)
- 树布局：Circular (ggtree)
- 图片分辨率：300 dpi
- 生成两个树版本：忽略分支长度 / 考虑分支长度

### 自动清理
脚本会自动删除旧版本的输出文件（以确保输出文件夹干净）

## 联系与反馈

如有问题或建议，请参考脚本注释或查看源代码。

---
*最后更新：2025-11-02*
