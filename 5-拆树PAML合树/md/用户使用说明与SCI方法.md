# mtDNA 拆树与合树 Pipeline 使用说明

## 1. 这个项目做什么

这个项目用于处理超过 `PAML baseml` 可高效处理规模的大型 mtDNA 系统发育树。

它完成两件事：

1. 拆树
   - 将大树拆分为多个适合 `baseml` 运行的 rooted 子树
   - 每棵 `baseml` 子树默认包含 `RSRS`
   - 小子树会自动补入局部锚点
2. 合树
   - 将每棵子树更新后的 branch lengths 回填到全局主树
   - 对非核心 scaffold 边进行启发式重标定

## 2. 目录结构

- `conf/Config.json`
  - 配置文件
- `data/`
  - 输入树和外群文件
- `pipe/run_split_pipeline.sh`
  - 拆树主控
- `pipe/run_merge_pipeline.sh`
  - 合树主控
- `python/`
  - Python 模块
- `output/`
  - 拆树输出
- `output/merge/`
  - 合树输出
- `md/`
  - 文档

## 3. 运行前准备

### 3.1 软件依赖

需要：

- `bash`
- `python3`
- `conda`
- `gotree`
- `Biopython`

当前项目默认要求：

- `gotree` 在 `conda` 环境中可用

环境检查脚本：
```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/script/check_env.sh
```

### 3.2 输入文件

必须准备：

- 原始树：`data/big.tree` 或 `data/big_rooted.tree`
- 外群 tip 文件：`data/outgroup_tips.txt`

当前工作流固定要求：

- 全局外群 tip 为 `RSRS(默认,也可在外群文件中指定)`
- 若输入树未定根，外群文件必须包含 `RSRS(默认,也可在外群文件中指定)`
- 若输入树已定根，则程序会直接使用已定根树

## 4. 配置文件说明

配置文件路径：[Config.json](/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/conf/Config.json)

关键参数如下。

### 4.1 拆树参数

- `paths.input_tree`
  - 输入树路径
- `paths.outgroup_tip_file`
  - 外群 tip 文件
- `paths.output_dir`
  - 输出目录
- `runtime.max_tips`
  - 单棵 `baseml` 子树最大 tip 数
- `runtime.min_baseml_tips`
  - 单棵 `baseml` 子树最小 tip 数
- `runtime.construct_baseml_subtrees`
  - 是否构建 `baseml` 子树
- `runtime.always_include_outgroup`
  - 是否始终加入 `RSRS`
- `runtime.reserve_slots_for_outgroup`
  - 从 `core max` 中预留给外群的位点数

### 4.2 合树参数

- `merge.enabled`
  - 是否启用合树流程
- `merge.analysis_tree_source`
  - `simulated` 或 `external`
- `merge.external_result_dir`
  - 真实 PAML 结果目录
- `merge.randomization_model`
  - 当前固定为 `lognormal_multiplicative`
- `merge.randomization_sigma`
  - 支长随机化强度
- `merge.randomization_seed`
  - 随机种子
- `merge.min_branch_length`
  - 最小正 branch length

## 5. 如何运行拆树

运行命令：
```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/run_split_pipeline.sh
```

### 5.1 拆树输出

主要输出位于 `output/`：

- `output/intermediate/rooted.tree`
  - 权威 rooted 主树
- `output/core_subtrees/`
  - 无重叠 `core` 子树
- `output/core_subtree_summary.tsv`
  - `core` 子树汇总
- `output/baseml_subtrees/`
  - 可直接用于 `baseml` 的子树
- `output/baseml_subtree_summary.tsv`
  - `baseml` 子树汇总
- `output/baseml_tree_manifest.tsv`
  - 每个 tip 在每棵子树中的角色
- `output/overlap_report.tsv`
  - 重叠统计
- `output/baseml_validation_report.tsv`
  - 拆树校验结果

### 5.2 如何理解两层子树

#### `core subtree`

特点：
- 全覆盖
- 无重叠
- 每个原始 tip 只属于一个 `core subtree`

用途：
- 全局映射
- 合树时确定唯一归属

#### `baseml subtree`

特点：

- 允许重叠
- 含 `RSRS`
- 可含 anchor tips

用途：

- 直接给 `PAML baseml`

## 6. 如何运行合树

### 6.1 模拟模式

如果还没有真实 PAML 结果，直接运行：

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/run_merge_pipeline.sh
```

当前默认会：
1. 对所有 `baseml` 子树做支长随机化
2. 模拟出一套“PAML 结果”
3. 回填到主树
4. 输出 merged tree

### 6.2 真实 PAML 结果模式

如果已经完成真实 `baseml` 计算：

1. 在 `Config.json` 中设置：

```json
"analysis_tree_source": "external"
```

2. 设置：

```json
"external_result_dir": "你的结果目录"
```

3. 结果目录中每棵树文件名必须是：

```text
baseml_subtree_0001.nwk
baseml_subtree_0002.nwk
...
```

4. 然后运行：

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/run_merge_pipeline.sh
```

### 6.3 合树输出

主要输出位于 `output/merge/`：

- `simulated_baseml_subtrees/`
  - 模拟结果树
- `simulation_manifest.tsv`
  - 模拟参数记录
- `merged_tree.nwk`
  - 最终合并树
- `merge_report.tsv`
  - 每个 core clade 的回填记录
- `edge_update_report.tsv`
  - 每条边的更新来源和尺度
- `merge_validation_report.tsv`
  - 合树校验报告
- `merge.log`
  - 合树日志

## 7. 合树结果怎么解释

最终 `merged_tree.nwk` 具有以下含义：

1. 拓扑与 rooted 主树保持一致
2. 对属于某个唯一 `core clade` 的边：
   - 使用对应子树结果中的 branch length
3. 对其余 scaffold 边：
   - 根据相邻已更新 core clades 的倍率做启发式重标定

因此 merged tree 不是简单地“拼接子树”，而是“主树模板 + core 边精确回填 + scaffold 尺度传播”。

## 8. 常见问题

### 8.1 输入树是否必须已定根

不是。

程序支持两种情况：
- 已定根：直接使用
- 未定根：使用外群文件重新定根

### 8.2 输入树是否允许多分叉

拆树递归本身允许内部多分叉，但当前完整 `baseml + merge` 工作流最终要求主树和结果树是 rooted binary tree。

### 8.3 为什么 `baseml subtree` 会重叠

因为小 `core subtree` 需要通过 anchor tips 和全局外群 `RSRS` 来构造可分析的 rooted 子树。重叠是设计允许的，不是错误。

### 8.4 为什么不能直接把整棵 `baseml subtree` 拼回去

因为 `baseml subtree` 中含有 anchor 和 `RSRS`，它们不是唯一属于该核心块的成员。直接整棵拼接会破坏全局唯一归属关系。

## 9. 推荐操作顺序

### 第一步：准备输入

确保：
- 输入树路径正确
- `RSRS` 在树中存在
- 外群文件正确

### 第二步：拆树

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/run_split_pipeline.sh
```

### 第三步：检查拆树结果

重点看：
- `output/baseml_validation_report.tsv`
- `output/baseml_subtree_summary.tsv`
- `output/overlap_report.tsv`

### 第四步：运行 PAML 或先做模拟

- 模拟：直接跑 merge pipeline
- 真实结果：先准备 `external_result_dir`

### 第五步：合树

```bash
bash /mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/pipe/run_merge_pipeline.sh
```

### 第六步：检查合树质量

重点看：

- `output/merge/merge_validation_report.tsv`
- `output/merge/merge_report.tsv`
- `output/merge/edge_update_report.tsv`

## 10. SCI 风格 Methods 示例

下面是一段可直接作为论文方法学草稿基础的中文版本。你后续可以根据期刊风格再精修。

### Methods

Large mitochondrial phylogenetic trees were partitioned into analysis-ready subtrees for downstream divergence-time inference using `PAML baseml`. The input phylogeny was provided in Newick format. If the input tree was already rooted, the original rooted topology was retained. Otherwise, the tree was rerooted using the outgroup sequence `RSRS` with `gotree reroot outgroup --strict`. A rooted master tree was then generated and used as the global topological reference for all downstream steps.

To enable large-tree analysis, the ingroup portion of the rooted master tree was recursively decomposed into non-overlapping monophyletic core partitions. Each core partition represented the largest monophyletic clade not exceeding the predefined size threshold. These core partitions were used as the unique ownership layer such that each original ingroup tip was assigned to exactly one core subtree. To construct `baseml`-ready analysis trees, each core partition was augmented with the global outgroup `RSRS`. When the resulting subtree size remained below the empirical minimum threshold for stable dating analysis, additional local anchor clades were iteratively added from phylogenetically adjacent sister lineages along the ancestral path. These anchor tips were allowed to overlap between analysis trees, whereas core partitions remained non-overlapping.

For each analysis subtree, the selected tip set was pruned from the rooted master tree to generate an induced rooted subtree. The induced subtree was then rerooted with `RSRS` and normalized to a strictly bifurcating topology. All generated `baseml` input trees were therefore rooted, binary, and explicitly contained the outgroup sequence.

To emulate branch-length updates prior to the availability of real `baseml` results, subtree branch lengths were randomly perturbed using a lognormal multiplicative model. Specifically, each non-root branch length was multiplied by a positive random factor sampled from a lognormal distribution with mean-preserving parameterization. These simulated trees were then treated as surrogate `baseml` outputs for testing the merge workflow.

For tree merging, the rooted master tree was used as the fixed global topology. Within each analysis tree, the clade corresponding to the original core partition was identified using the exact set of core tips. Branch lengths for all edges uniquely belonging to that core clade, including terminal edges, internal edges, and the clade entry edge, were mapped back to the corresponding edges in the master tree using descendant tip-set signatures. For scaffold edges that were not uniquely covered by a single core clade, branch lengths were not copied directly from any single subtree. Instead, a heuristic rescaling factor was estimated from adjacent updated core clades. The scaling factor for each scaffold edge was computed as the `core_n_tips`-weighted geometric mean of descendant core-root multipliers, and the original scaffold branch length was rescaled accordingly. This produced a merged tree that preserved the original global topology while incorporating subtree-specific branch-length updates and globally coherent branch-length scaling.

The final merged tree was validated against the rooted master tree to ensure identical tip membership, identical tip naming, preservation of rooted bifurcating topology, and correct placement of `RSRS` as the singleton outgroup child of the root. Only branch lengths were allowed to differ between the merged tree and the rooted master tree.

## 11. 建议引用的软件

建议在 Methods 或 Software 段落中明确写出：

- `PAML baseml`
- `gotree`
- `Biopython`

