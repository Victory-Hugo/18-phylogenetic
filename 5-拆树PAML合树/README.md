# mtDNA Backbone-Target-Graft Pipeline

## 1. 流程概述

本项目当前采用 `骨架树 + 目标 clade graft` 。

完整流程为：

```sh
输入大树 + FASTA + 外群列表
-> rooted 主树确认
-> 确定 backbone tip 集
-> 生成 backbone tree
-> 递归切分 target subtree
-> 构造 PAML 子树 = 全 backbone + target + outgroup
-> 准备 treefile / fasta / ctl
-> 批量运行 baseml
-> 解析 analysis_trees
-> 构建 assembly_scaffold
-> 聚合 backbone 边长
-> graft 各 target clade
-> 输出 merged_tree.nwk
```

这套逻辑的关键原则是：

- 每一棵 PAML 子树都保留同一套骨架树
- 目标 clade 只负责自己新增的非骨架 tips

## 2. 当前配置默认值

分阶段配置文件：

- [1-split.yaml](/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/conf/1-split.yaml)
- [2-paml.yaml](/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/conf/2-paml.yaml)
- [3-merge.yaml](/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/conf/3-merge.yaml)
- [4-ultrastandard.yaml](/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/conf/4-ultrastandard.yaml)

当前默认样例输入：

- 树：`data/big_rooted.tree`
- FASTA：`data/big.fasta`
- 外群：`conf/outgroup_tips.txt`
- `runtime.max_tips = 300`
- `runtime.backbone_size = 100`
- `runtime.backbone_sampling_strategy = greedy_farthest_patristic`
- `runtime.target_partition_mode = recursive_monophyletic`
- `paml.seq_id_strategy = prefix_before_first_underscore`
- `merge.mode = backbone_graft`

## 3. 输入要求

### 3.1 树文件

- 格式为 `Newick`
- 可以是已定根树
- 若未定根，必须提供外群列表，由 `gotree reroot outgroup --strict` 重新定根

### 3.2 FASTA 文件

- 格式为 `FASTA`
- 必须能与树 tip 对应


### 3.3 外群文件

- UTF-8 纯文本
- 每行一个 tip
- 支持空行和 `#` 注释

当前样例必须包含：

```text
RSRS_MRCA
```

若外群文件中只有这 1 个 tip，则无需再在 YAML 中重复设置 `outgroup_tip_name`。

## 4. 自动骨架树逻辑

骨架来源优先级：

1. `paths.backbone_tip_id_file`
2. `paths.backbone_tree`
3. 自动抽样

自动抽样不是简单随机抽样，而是：

1. 先把 ingroup 递归切成 `backbone_size` 个前沿 clade
2. 再从每个前沿 clade 中选择一个代表 tip
3. 代表 tip 默认选取该前沿 clade 中距离前沿根最深的终端 tip

这样做的目的，是避免 backbone tip 过度分散到太多局部小分支，减少 target clade 被碎片化。

输出文件：

- `output/split/backbone_tips.txt`
- `output/split/backbone_tree.nwk`
- `output/split/backbone_summary.tsv`

## 5. target subtree 逻辑

当前 target subtree 不再等同于旧版 core subtree。

定义规则：

- 只在 rooted master tree 的 ingroup 部分递归切分
- 切分依据是“非骨架 tip 数”
- 一个 clade 只要其非骨架 tip 数 `<= target_capacity`，即可直接作为 target subtree
- 该 clade 可以含有 backbone descendant；真正用于 graft 的对象是该 clade 内目标 tip 的诱导子树，而不是把整块 backbone 一起替换掉

其中：

```text
target_capacity = max_tips - backbone_size - 1
```

最后的 `1` 用于给全局外群保留位置。

输出文件：

- `output/split/target_subtree_summary.tsv`
- `output/split/target_tree_manifest.tsv`
- `output/split/target_subtrees/*.nwk`

说明：

- `target_subtree_summary.tsv` 记录的是该 target 负责的非骨架 tips
- `target_subtrees/*.nwk` 保存的是对应 target clade 的原始子树，因此可能同时包含少量 backbone descendant

## 6. PAML 子树逻辑

每棵 PAML 子树的 tip 集固定为：

```text
backbone tips + target non-backbone tips + outgroup tip
```

也就是说，骨架树始终完整保留，target 只额外带自己负责的那一部分样本。

输出文件：

- `output/split/paml_subtree_summary.tsv`
- `output/split/paml_tree_manifest.tsv`
- `output/split/paml_subtrees/*.nwk`

## 7. PAML 输入准备规则

### 7.1 treefile

每棵 `treefile` 强制写成两行：

```text
  <tip数> 1
<单行 Newick>
```

### 7.2 FASTA

每棵子树 FASTA 的 header 必须与对应树 tip 完全一致。

### 7.3 ctl

`ctl` 采用模板改写方式：

- 只重写 `seqfile`
- 只重写 `treefile`
- 只重写 `outfile`

其他参数保持模板原样。

当前实现不是在 Python 中硬编码整份 `baseml.ctl` 内容，而是从 `conf/2-paml.yaml` 的 `paml.ctl_template` 读取模板；默认模板文件为 `conf/baseml_mtDNA.ctl`。

因此：

- 若要修改 `model`、`clock`、`kappa`、`alpha`、`cleandata` 等 PAML 参数，直接修改模板文件即可
- 程序只会为每个子树 job 改写 3 个字段：`seqfile`、`treefile`、`outfile`
- 若模板中缺少这 3 个字段，程序会在生成的 `baseml.ctl` 末尾自动补上

## 8. merge 逻辑

merge 的核心是 scaffold + graft：

1. 从原始 rooted master tree 生成 `assembly_scaffold.nwk`
2. scaffold 只保留：
   - backbone tips
   - target placeholders
   - outgroup
3. 从每棵 analysis tree 中提取对应 target clade
4. 先聚合 backbone 边长
5. 再把 target clade graft 回 scaffold

输出文件：

- `output/merge/assembly_scaffold.nwk`
- `output/merge/backbone_edge_estimates.tsv`
- `output/merge/graft_report.tsv`
- `output/merge/edge_update_report.tsv`
- `output/merge/merged_tree.nwk`
- `output/merge/merge_validation_report.tsv`

## 9. 运行方式

四个阶段完全独立，分别读取自己的 YAML 配置。

拆树：

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/1-run_split_pipeline.sh
```

PAML：

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/2-run_paml_pipeline.sh
```

合树：

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/3-run_merge_pipeline.sh
```

超度量树标准化：

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/4-run_ultrastandard_pipeline.sh
```

如需切换配置，可显式指定：

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/1-run_split_pipeline.sh --config /path/to/1-split.yaml
```

## 10. 注意

骨架树的选择对结果影响较大，**建议优先使用** `paths.backbone_tip_id_file` 明确指定骨架，而不是完全依赖自动抽样。