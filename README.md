# 18-phylogenetic

![phylo.png|550](https://picturerealm.oss-cn-chengdu.aliyuncs.com/obsidian/20260318181321893.png)

mtDNA 系统发育分析工具与流程集合，覆盖 PAML 分子钟估算、Rho 定年、pplacer 样本置放、BEAST2 贝叶斯分析，以及“大树拆分 -> 批量 PAML -> 合并超度量树”的主流程。

## 当前目录结构

```text
18-phylogenetic/
├── 0-PAML_v10_Linux/   PAML v10 Linux 预编译包
├── 1-PAML结果整理/      小中型数据的标准 baseml 流程
├── 2-Rho/              Rho 方法估算 MRCA 时间
├── 3-pplacer/          参考树置放与结果可视化
├── 4-BEAST2/           BEAST2 软件包与示例
├── 5-拆树PAML合树/     大树分治 + PAML + 合并 + 超度量/时间校正主流程
└── 6-骨架树筛选/        为主流程筛选多层级 backbone 样本
```

## 模块说明

### 0-PAML_v10_Linux

PAML v10 的 Linux 预编译目录，包含常规 `baseml`，以及一个自定义的多线程版本 `baseml_llt`。

- 适用场景：本地直接运行 PAML，或给 `5-拆树PAML合树` 提供可执行程序
- 注意：`baseml_llt` 为自定义版本，不是官方发布版，建议用于测试和性能评估

示例：

```bash
BASEML_PAR3=1 OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores \
  0-PAML_v10_Linux/bin/baseml_llt your_baseml.ctl
```

### 1-PAML结果整理

标准 `baseml` 批处理目录，适合中小规模树或单次结果整理。

当前 `pipe/` 中的主要脚本：

- `1-tree-number统计.sh`：统计输入树 tip 数量
- `2_运行.sh`：批量运行 `baseml`
- `3_结果整理.sh`：汇总 PAML 输出

目录自带 `data/`、`meta/`、`examples/` 和 `conf/`，可作为最小 baseml 示例工程使用。

### 2-Rho

使用 Rho 方法估算单倍群最近共同祖先时间。

当前流程入口位于 `pipe/`：

- `1_准备单倍群列表.ipynb`
- `2-准备序列文件.sh`
- `3-vcf→fasta_通过plink.sh`
- `4_RHO方法核心.sh`
- `5_RHO计算时间.sh`
- `9_结果整理.ipynb`

详细说明见 [2-Rho/README.md](2-Rho/README.md)。

### 3-pplacer

将待测序列置放到参考系统发育树，适合大规模样本的快速定位。

当前脚本位于 `script/`：

- `1-建树→pplacer.sh`：构建参考树并运行 `pplacer`
- `2-1-可视化.sh`：可视化辅助脚本
- `2-可视化.R`：绘制置放结果

补充说明见 [3-pplacer/script/README_可视化.md](3-pplacer/script/README_可视化.md)。

### 4-BEAST2

BEAST2 软件目录，包含程序本体、JRE、示例 XML 和启动资源。

- 适用场景：贝叶斯分子钟、分歧时间估算、区间推断
- 当前目录更偏向软件分发与示例保留，不是统一封装后的 shell 流程

### 5-拆树PAML合树

当前主流程。目标是在 PAML 无法直接处理超大树时，通过固定骨架树把大树拆分成多个子树，批量运行 `baseml` 后再合并回全局树，并继续做超度量化与时间校正。

主流程当前是 5 个阶段：

1. `split`：拆树并设计 PAML 子树
2. `paml`：准备输入、批量运行并解析输出
3. `merge`：把各子树结果 graft 回主树
4. `ultrastandard`：投影为超度量树
5. `time calibration`：将支长换算为绝对时间

一键全流程：

```bash
bash 5-拆树PAML合树/pipe/0-run_all_pipeline.sh
```

分阶段运行：

```bash
bash 5-拆树PAML合树/pipe/1-run_split_pipeline.sh
bash 5-拆树PAML合树/pipe/2-run_paml_pipeline.sh
bash 5-拆树PAML合树/pipe/3-run_merge_pipeline.sh
bash 5-拆树PAML合树/pipe/4-run_ultrastandard_pipeline.sh
bash 5-拆树PAML合树/pipe/5-run_time_calib_pipeline.sh
```

当前配置文件：

- `conf/1-split.yaml`
- `conf/2-paml.yaml`
- `conf/3-merge.yaml`
- `conf/4-ultrastandard.yaml`
- `conf/5-time_calib.yaml`

详细说明见 [5-拆树PAML合树/README.md](5-拆树PAML合树/README.md)。

### 6-骨架树筛选

为 `5-拆树PAML合树` 提供 backbone 候选样本筛选。当前实现基于单倍群层级覆盖和变异距离，输出多套嵌套的 backbone 名单。

当前入口：

```bash
bash 6-骨架树筛选/pipe/run_backbone_selection.sh
```

默认配置文件：

- `conf/Config.yaml`

主要输入：

- `meta/骨架树筛选.tsv`
- `meta/基础信息.tsv`
- `data/merged_clean_去除hot.snp.vcf.gz`

主要输出：

- `output/backbone_100.tsv`
- `output/backbone_150.tsv`
- `output/backbone_200.tsv`
- `output/backbone_250.tsv`
- `output/backbone_300.tsv`
- `output/05_backbone_selection_master.tsv`
- `output/README_筛选说明.md`

## 推荐使用顺序

如果目标是跑“大树拆分 + PAML 合并”主方案，建议按下面顺序使用：

1. 先运行 `6-骨架树筛选`，确定 backbone 样本集合
2. 将 backbone tip 列表写入 `5-拆树PAML合树/conf/` 或对应输入路径
3. 运行 `5-拆树PAML合树/pipe/0-run_all_pipeline.sh`
4. 按需使用 `2-Rho`、`3-pplacer` 或 `4-BEAST2` 做交叉验证或补充分析

## 常用依赖

| 工具 | 用途 |
|------|------|
| [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) | 最大似然分子钟估算 |
| [BEAST2](https://www.beast2.org/) | 贝叶斯分子钟分析 |
| [pplacer](https://matsen.fhcrc.org/pplacer/) | 样本置放到参考树 |
| [gotree](https://github.com/evolbioinfo/gotree) | 树定根、检查、剪枝与格式处理 |
| Python >= 3.10 | 主流程与辅助脚本 |
| R | pplacer 可视化 |

## 备注

- 当前仓库里既有“软件包目录”（如 `0-PAML_v10_Linux`、`4-BEAST2`），也有“可直接运行的流程目录”（如 `2-Rho`、`5-拆树PAML合树`、`6-骨架树筛选`）。
- 顶层 README 只做导航；具体参数、输入约束和输出解释以各子目录 README 或 `conf/*.yaml` 为准。
