# Backbone-Target-Graft 代码逻辑说明

## 1. 设计目标

当前代码已经从旧版 `core subtree` 方案切换为 `backbone tree + target graft` 方案。新的目标不是把大树机械切成很多互不相交的小块，而是：

- 让每棵 PAML 儿子树都保留足够大的全局遗传多样性背景
- 让每棵儿子树只额外负责一个局部 target clade
- 让最终合树通过 scaffold + graft 完成，而不是依赖旧版 core 回填

## 2. 主控脚本

### `pipe/1-run_split_pipeline.sh`

职责：

1. 读取 `conf/1-split.yaml`
2. 解析输入树、外群、backbone 参数
3. 调用 [split_phylo_tree.py](/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/python/split_phylo_tree.py)
4. 调用 [validate_phylo_split.py](/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/python/validate_phylo_split.py)

### `pipe/2-run_paml_pipeline.sh`

职责：

1. 读取 `conf/2-paml.yaml`
2. 读取 `output/split/paml_subtree_summary.tsv`
2. 生成 `jobs/*/treefile, fasta, ctl`
3. 批量运行 baseml
4. 解析 `mlb`

### `pipe/3-run_merge_pipeline.sh`

职责：

1. 读取 `conf/3-merge.yaml`
2. 读取新的 split 输出
3. 加载 analysis trees
4. 构建 scaffold
5. 聚合 backbone 边长
6. graft 各 target clade
7. 校验 merged tree

## 3. split 逻辑

### 3.1 rooted tree 确认

公共函数在 [phylo_split_common.py](/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/python/phylo_split_common.py)：

- `detect_tree_rooted_status`
- `prepare_rooted_tree`
- `get_root_children_for_outgroup`

要求：

- root 下必须存在唯一 singleton outgroup child
- 若 `outgroup_tip_file` 仅包含 1 个 tip，则可自动推断 singleton outgroup，无需额外设置 `outgroup_tip_name`

### 3.2 backbone 抽样

自动 backbone 抽样在 [backbone_sampling.py](/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/python/backbone_sampling.py)。

当前实现不是把 100 个 tip 直接用全局 tip 级 farthest-point 撒到整棵树上，而是：

1. 从 ingroup root 开始
2. 递归分裂最大 clade，直到得到 `backbone_size` 个前沿 clade
3. 每个前沿 clade 只选 1 个代表 tip
4. 代表 tip 默认取该前沿 clade 内最深的终端 tip

这样可以减少骨架 tip 在局部区域过度堆积的问题，降低 target clade 过碎的风险。

### 3.3 target subtree 切分

函数：

- `compute_target_partition_counts`
- `build_target_partition`

规则：

- 只在 ingroup 部分递归
- 容量约束按“非骨架 tip 数”计算
- 允许 target 所在 clade 含有 backbone descendant
- scaffold 与 graft 始终按 target tip 的诱导子树操作，而不是整块原 clade 直接替换

公式：

```text
target_capacity = max_tips - backbone_size - 1
```

其中 `1` 是全局 outgroup 的位置。

### 3.4 PAML subtree 构造

函数：

- `build_paml_subtrees`
- `build_induced_tree_file`

每棵子树 tip 集固定为：

```text
backbone tips + target non-backbone tips + outgroup tip
```

输出文件：

- `backbone_summary.tsv`
- `target_subtree_summary.tsv`
- `target_tree_manifest.tsv`
- `paml_subtree_summary.tsv`
- `paml_tree_manifest.tsv`

## 4. PAML 层逻辑

### `prepare_paml_inputs.py`

职责：

- 从 `paml_subtree_summary.tsv` 读取每棵子树的 tip 集
- 根据 `seq_id_strategy` 从总 FASTA 中提取对应序列
- 重写 FASTA header，使其与树 tip 完全一致
- 生成两行格式的 `treefile`
- 渲染 `baseml.ctl`

### `run_paml_jobs.py`

职责：

- 批量运行 baseml
- 支持并行
- 支持 `skip_existing`
- 支持 `fh negative` 时 fallback 到 `fix_blength = 2`

### `parse_paml_outputs.py`

职责：

- 解析 `mlb`
- 提取 named tree
- 生成 `analysis_trees/*.nwk`
- 生成 `Parameters.csv`
- 生成 `TIP_Length.csv`

## 5. merge 逻辑

merge 公共函数位于 [phylo_merge_common.py](/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/python/phylo_merge_common.py)。

### 5.1 scaffold 构建

核心函数：

- `build_scaffold_tree`

逻辑：

1. 从原始 rooted master tree 出发
2. 对每个 target subtree root 用 `TARGET_xxxx` 替换
3. 保留 backbone tips 和 outgroup
4. 生成 `assembly_scaffold.nwk`

### 5.2 backbone 边长聚合

核心函数：

- `compute_branch_by_signature`
- `build_backbone_signature_key`
- `weighted_geometric_mean`

逻辑：

1. 从每棵 analysis tree 中提取 backbone 对应边
2. 按 descendant backbone tip signature 匹配
3. 对同一 signature 的边长度做加权几何平均
4. 回写到 scaffold 主干

### 5.3 target graft

核心函数：

- `extract_monophyletic_target_clade`
- `replace_placeholder_with_clade`

逻辑：

1. 用 `target_nonbackbone_tip_names` 在 analysis tree 中找 MRCA
2. 要求该 target tip 集严格单系
3. 用提取出的 clade 替换 scaffold 中对应的 `TARGET_xxxx`

### 5.4 merged tree 校验

在 [validate_merged_tree.py](/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/python/validate_merged_tree.py) 中检查：

- scaffold 只含 backbone tips、placeholders、outgroup
- merged tree 不残留任何 placeholder
- merged tree tip 集与 master tree 一致
- rooted
- binary
- 所有 graft 均成功记录

## 6. 旧版逻辑已废弃

以下旧概念不应再出现在主流程代码推理中：

- `core_subtree`
- `anchor_tip`
- `min_baseml_tips`
- `strict_core_topology_match`
- `scaffold_scale_aggregation`

如果未来维护时又看到这些字段，默认视为历史遗留，不应继续作为当前设计前提。

## 7. 当前实现的现实约束

当前版本解决了你指出的关键问题：PAML 子树不再缺少全局遗传多样性背景。

但仍有一个现实 tradeoff：

- 如果 backbone tip 过于分散，target subtree 数量仍然会上升
- 为了缓解这一点，自动 backbone 抽样已经从“全局散点”改为“前沿 clade 代表 tip”策略

这意味着当前实现的优化方向已经明确：

1. 先减少 backbone 对局部 clade 的碎裂作用
2. 再在 target partition 上尽量保持单系与可 graft 性

如果后续还要继续优化，优先改的是：

- backbone 代表 tip 的选择策略
- target partition 的后处理聚合策略

而不是回退到旧版 `core + anchor` 方案。
