# 18-phylogenetic

mtDNA 系统发育分析流程集合，涵盖系统发育树构建、分子钟估算与定年等核心任务。

## 目录结构

```
18-phylogenetic/
├── 0-PAML_v10_Linux/   PAML v10 Linux 二进制包（含自定义多线程版 baseml_llt）
├── 1-PAML/             标准 PAML baseml 分析流程
├── 2-Rho/              Rho 方法估算 MRCA 时间
├── 3-pplacer/          pplacer 系统发育样本置放
├── 4-BEAST2/           BEAST2 贝叶斯分子钟分析
└── 5-拆树PAML合树/     大树拆分→批量 PAML→合并超度量树 主流程
```

---

## 模块说明

### 0-PAML_v10_Linux

PAML v10 的 Linux 预编译包，同时包含为 `baseml` 定制的 **OpenMP 多线程版本** `baseml_llt`。

>请注意，这不是一个官方发布的 PAML 版本，且其功能不稳定，请勿在生产环境中使用。建议仅用于测试和性能评估。

- 详细编译与运行说明见 [0-PAML_v10_Linux/BASEML_LLT.md](0-PAML_v10_Linux/BASEML_LLT.md)

推荐运行方式：

```bash
BASEML_PAR3=1 OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores \
    0-PAML_v10_Linux/bin/baseml_llt your_baseml.ctl
```

---

### 1-PAML

标准 `baseml` 分析流程，适用于中小规模数据集的直接运行。

流程脚本（`pipe/`）：

| 脚本 | 功能 |
|------|------|
| `1-tree-number统计.sh` | 统计输入树的 tip 数量 |
| `2_运行.sh` | 批量运行 baseml |
| `3_结果整理.sh` | 汇总 baseml 输出结果 |

详细说明见 [1-PAML/README.md](1-PAML/README.md)。

---

### 2-Rho

使用 **Rho 方法**（ρ 法）估算 mtDNA 单倍群的最近共同祖先（MRCA）时间。

流程脚本（`pipe/`）：

| 脚本 | 功能 |
|------|------|
| `1_准备单倍群列表.ipynb` | 生成目标单倍群列表 |
| `2-准备序列文件.sh` | 从 VCF/FASTA 中提取子集 |
| `3-vcf→fasta_通过plink.sh` | VCF 转 FASTA（plink 辅助） |
| `4_RHO方法核心.sh` | 运行 Rho 核心计算 |
| `5_RHO计算时间.sh` | 将 ρ 值换算为绝对时间 |
| `9_结果整理.ipynb` | 整理和可视化结果 |

详细说明见 [2-Rho/README.md](2-Rho/README.md)。

---

### 3-pplacer

使用 **pplacer** 将待测序列置放到参考系统发育树上，适用于大规模样本的快速置放场景。

脚本（`script/`）：

| 脚本 | 功能 |
|------|------|
| `1-建树→pplacer.sh` | 构建参考树并运行 pplacer |
| `2-可视化.R` | R 可视化置放结果 |
| `2-1-可视化.sh` | Shell 辅助可视化 |

---

### 4-BEAST2

**BEAST2** 贝叶斯分子钟分析，用于估算系统发育树的分歧时间与置信区间。

---

### 5-拆树PAML合树（主流程）

核心流程，解决 PAML 无法直接处理超大树的问题。将大树拆分为若干含共享骨架的子树，批量运行 PAML，再将结果合并还原为一棵超度量树。

**流程概览：**

```
输入大树 + FASTA + 外群列表
  → 确定骨架树（backbone）
  → 递归切分 target subtree
  → 构造 PAML 子树（全骨架 + target + 外群）
  → 批量运行 baseml
  → 聚合骨架边长 + graft target clade
  → 输出 merged_tree.nwk
```

**四阶段独立运行：**

```bash
# 1. 拆树
bash 5-拆树PAML合树/pipe/1-run_split_pipeline.sh

# 2. 批量 PAML
bash 5-拆树PAML合树/pipe/2-run_paml_pipeline.sh

# 3. 合树
bash 5-拆树PAML合树/pipe/3-run_merge_pipeline.sh

# 4. 超度量树标准化
bash 5-拆树PAML合树/pipe/4-run_ultrastandard_pipeline.sh
```

**配置文件（`conf/`）：**

| 文件 | 对应阶段 |
|------|----------|
| `1-split.yaml` | 拆树参数（骨架大小、target 分区策略等） |
| `2-paml.yaml` | PAML 运行参数与 ctl 模板路径 |
| `3-merge.yaml` | 合树模式与 graft 参数 |
| `4-ultrastandard.yaml` | 超度量树标准化参数 |

详细说明见 [5-拆树PAML合树/README.md](5-拆树PAML合树/README.md)。

---

## 依赖工具

| 工具 | 用途 |
|------|------|
| [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) | 最大似然分子钟估算 |
| [BEAST2](https://www.beast2.org/) | 贝叶斯分子钟分析 |
| [pplacer](https://matsen.fhcrc.org/pplacer/) | 系统发育样本置放 |
| [gotree](https://github.com/evolbioinfo/gotree) | 树操作（定根、剪枝等） |
| Python ≥ 3.10 | 流程脚本 |
| R | pplacer 结果可视化 |
