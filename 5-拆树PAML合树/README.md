# mtDNA Backbone-Target-Graft Pipeline

> 说明：本目录是从 `7-PAML对ML进行校正` 复制出的“新方案开发基线工程”，用于后续优化“PAML 大树分治后再合并”的方法学流程。
>
> 本轮仅完成工程基线复制与长期记忆文档落地，不实现新算法。`output/` 未被复制，后续需要在本目录内重新运行各阶段流水线产生结果。

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
-> 超度量标准化
-> 时间校正（分子钟 / 节点标记）
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
- [5-time_calib.yaml](/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/conf/5-time_calib.yaml)

当前默认样例输入：

- 树：`data/big_rooted.tree`
- FASTA：`data/big.fasta`
- 外群：`conf/outgroup_tips.txt`
- `runtime.max_tips = 300`
- `runtime.backbone_size = 100`
- `runtime.backbone_sampling_strategy = greedy_farthest_patristic`
- `runtime.target_partition_mode = recursive_monophyletic`
- `paml.seq_id_strategy = prefix_before_first_underscore`（优先完整匹配，找不到时再回退到首个下划线前缀）
- `merge.mode = backbone_graft`
- `ultrastandard.projection_mode = constrained_height_fit`
- `time_calib.method = molecular_clock`
- `time_calib.substitution_rate_per_site_per_year = 2.67e-8`
- `time_calib.sequence_length = 16569`

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
RSRS
```

若外群文件中只有这 1 个 tip，则无需再在 YAML 中重复设置 `outgroup_tip_name`。

## 4. 自动骨架树逻辑

骨架来源优先级：

1. `paths.backbone_tip_id_file`
2. 自动抽样

兼容旧配置时，脚本仍接受 `paths.backbone_tree`，但默认 `conf/1-split.yaml` 已不再保留这个字段。

当提供 `paths.backbone_tip_id_file` 时，当前行为是：

- 若文件中的 tip 数量恰好等于 `runtime.backbone_size`，则直接使用
- 若文件中的 tip 数量小于 `runtime.backbone_size`，则保留这些用户指定 tip，并继续自动抽样补齐到 `runtime.backbone_size`
- 若文件中的 tip 数量大于 `runtime.backbone_size`，则流程会先发出警告，再终止运行

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
target_capacity = max_tips - 实际 backbone tip 数 - local_anchor_count - 1
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
backbone tips + target non-backbone tips + local anchors + outgroup tip
```

也就是说，骨架树始终完整保留，target 只额外带自己负责的那一部分样本，并保留局部 overlap anchors 作为后续全局尺度统一的共享参照。

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

## 9. ultrastandard 逻辑

ultrastandard 阶段在不改变树拓扑的前提下，将 `merged_ml_tree.nwk` 投影为超度量树。

当前支持两种模式：

1. `extend_terminal_to_max_depth`
   - 仅延长终端枝到最大 root-to-tip 深度
   - 内部枝长尽量保持不变
2. `constrained_height_fit`
   - 在最小正枝长约束下重估节点高度
   - 可同时调整终端枝和内部枝

输出文件：

- `output/ultrastandard/projection_input_tree.nwk`
- `output/ultrastandard/ultrametric_projection_report.tsv`
- `output/ultrastandard/root_to_tip_report.tsv`
- `output/ultrastandard/merged_ultrametric_tree.nwk`
- `output/ultrastandard/ultrastandard_validation_report.tsv`

## 10. 时间校正逻辑

时间校正阶段默认读取：

- `output/ultrastandard/merged_ultrametric_tree.nwk`

并将支长整体缩放到“年”为单位。

当前支持两种方法：

### 10.1 分子钟法

配置：

- `time_calib.method = molecular_clock`

默认参数：

- `time_calib.substitution_rate_per_site_per_year = 2.67e-8`
- `time_calib.sequence_length = 16569`
- `time_calib.branch_length_unit = substitutions_per_site`

说明：

- 若支长单位是“每位点替换数”，则按 `branch_length / rate` 换算为年
- 若将来输入树支长单位是“整条序列总替换数”，可改为 `substitutions_per_sequence`，此时会额外使用 `sequence_length` 参与换算

### 10.2 节点标记法

配置：

- `time_calib.method = node_calibration`

默认参数：

- `time_calib.node_calibration_tip_name = RSRS`
- `time_calib.node_calibration_divergence_years = 180000`

说明：

- 该方法会把 `RSRS` 这条谱系相对其余样本的分歧时间固定为距今 18 万年
- 然后据此对整棵超度量树做整体缩放

输出文件：

- `output/time_calib/merged_ultrametric_tree_years.nwk`
- `output/time_calib/time_calibration_summary.tsv`
- `output/time_calib/time_calibration_edge_report.tsv`

## 11. 运行方式

五个阶段完全独立，分别读取自己的 YAML 配置。

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

时间校正：

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/5-run_time_calib_pipeline.sh
```

一键串行运行 1-5 步：

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/0-run_all_pipeline.sh
```

如需切换配置，可显式指定：

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/1-run_split_pipeline.sh --config /path/to/1-split.yaml
```

总控脚本也支持分别指定 5 个配置文件：

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/0-run_all_pipeline.sh \
  --split-config /path/to/1-split.yaml \
  --paml-config /path/to/2-paml.yaml \
  --merge-config /path/to/3-merge.yaml \
  --ultrastandard-config /path/to/4-ultrastandard.yaml \
  --time-calib-config /path/to/5-time_calib.yaml
```

## 12. 注意

骨架树的选择对结果影响较大，**建议优先使用** `paths.backbone_tip_id_file` 明确指定骨架，而不是完全依赖自动抽样。

时间校正阶段只是对超度量树做整体尺度换算，不会改变树拓扑。若年代明显偏大或偏小，优先检查：

- 上游 `merged_ultrametric_tree.nwk` 的 root-to-tip 深度是否合理
- 分子钟速率是否适用于当前数据
- 是否应改用 `node_calibration` 而不是 `molecular_clock`
