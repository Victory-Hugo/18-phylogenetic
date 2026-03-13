# mtDNA 拆树与合树 Pipeline 代码逻辑说明

## 1. 项目目标
本项目服务于大规模 mtDNA 系统发育树的 `PAML baseml` 分块定年分析。整体流程分为两部分：

1. 拆树：把大树拆成适合 `baseml` 运行的多个 rooted 子树。
2. 合树：在得到各子树的定年结果后，把更新后的 branch lengths 回填到全局主树，并对未被直接更新的 scaffold 边做启发式重标定。

本说明面向人工智能或后续维护代理，重点解释代码职责、模块依赖、关键数据结构和控制流。

## 2. 目录职责
- `conf/Config.json`
  - 唯一配置入口。
- `pipe/run_split_pipeline.sh`
  - 拆树主控。
- `pipe/run_merge_pipeline.sh`
  - 合树主控。
- `script/check_env.sh`
  - 环境检查。
- `python/config_loader.py`
  - 从 `Config.json` 读取单个键值。
- `python/phylo_split_common.py`
  - 拆树与 `baseml` 子树构造的公共函数库。
- `python/split_phylo_tree.py`
  - 拆树主模块。
- `python/validate_phylo_split.py`
  - 拆树结果校验模块。
- `python/phylo_merge_common.py`
  - 合树公共函数库。
- `python/simulate_baseml_results.py`
  - 模拟 PAML 结果。
- `python/merge_baseml_subtrees.py`
  - 合并子树结果回主树。
- `python/validate_merged_tree.py`
  - 合树结果校验。
- `output/`
  - 拆树阶段产物。
- `output/merge/`
  - 合树阶段产物。

## 3. 配置读取逻辑
### `python/config_loader.py`
职责非常单一：

1. 读取 `Config.json`
2. 解析点号路径键，例如 `paths.input_tree`
3. 输出对应值

这个模块被两个 shell 主控反复调用，因此它不承担业务逻辑，只承担“配置桥接”功能。

关键约束：
- `run(config, key)` 返回字符串
- `main(argv=None)` 负责 CLI
- 不读取任何隐式运行状态

## 4. 拆树模块逻辑
### 4.1 拆树主控：`pipe/run_split_pipeline.sh`
控制流如下：

1. 定位项目根目录。
2. 调用 `script/check_env.sh`。
3. 使用 `config_loader.py` 逐项读取配置。
4. 创建输出目录。
5. 如果 `clean_subtree_dir=true`，清除旧结果。
6. 调用 `python/split_phylo_tree.py`。
7. 调用 `python/validate_phylo_split.py`。

这个 shell 只做调度，不做树算法。

### 4.2 公共拆树库：`python/phylo_split_common.py`
这个文件是拆树逻辑的核心。

#### A. 基础 I/O 与日志
- `setup_logger()`
- `read_newick_tree()`
- `write_tree()`
- `write_table()`
- `load_table()`
- `write_root_list()`

#### B. 输入预处理
- `parse_tip_file()`
- `validate_outgroup_tips()`
- `load_backbone_tips()`
- `detect_tree_rooted_status()`
- `prepare_rooted_tree()`

这里的关键点是：
- 若输入树已经 rooted，直接复制为 `output/intermediate/rooted.tree`
- 若输入树 unrooted，则必须用外群文件经 `gotree reroot outgroup --strict` 定根

#### C. 树结构索引
- `assign_node_ids()`
  - 给每个内部节点分配稳定 `NODE_XXXXXX`
  - 给 tip 分配 `TIP::<tip_name>`
- `compute_tip_counts()`
  - 计算每个 clade 的后代 tip 数

后续的拆树、校验、合树都依赖这套编号。

#### D. core partition 构造
- `collect_selected_roots_from_nodes()`
- `merge_small_blocks()`
- `build_core_partition()`

逻辑是：

1. 先确认根的两个子支中，一个必须是仅包含 `RSRS` 的 singleton outgroup。
2. 只对 ingroup 部分做核心分区。
3. 用“最大但不超过阈值的单系块递归拆分”选出 `core clade`。
4. 每个原始 ingroup tip 必须且只能属于一个 `core subtree`。

#### E. `baseml` 子树构造
- `gather_candidate_anchor_units()`
- `build_baseml_subtrees()`
- `build_induced_tree_file()`

逻辑是：

1. 每个 `core subtree` 先加上全局外群 `RSRS`。
2. 如果总 tip 数仍小于 `min_baseml_tips`，则沿祖先路径搜索姐妹 clade 作为锚点。
3. 锚点允许重叠，但 `core` 层不允许重叠。
4. 最终的 `baseml subtree` 不是某个单一 clade，而是从 rooted master tree 中按选定 tip 集修剪得到的 induced subtree。

#### F. 二叉化与合法性
- `normalize_tree_binary()`
- `is_binary_rooted_with_outgroup()`

这一步强制 `baseml subtree` 成为 rooted binary tree，并保证 `RSRS` 是根下 singleton child。

### 4.3 拆树入口：`python/split_phylo_tree.py`
这个模块是拆树主逻辑的组装层。

执行顺序：

1. 读取输入树。
2. 检查 `RSRS` 是否存在。
3. 判断树是否已定根。
4. 生成 `output/intermediate/rooted.tree`。
5. 对 rooted tree 重新编号。
6. 构建 `core partition`。
7. 输出：
   - `core_subtree_summary.tsv`
   - `core_subtree_roots.txt`
   - `core_subtrees/*.nwk`
8. 构建 `baseml subtrees`。
9. 输出：
   - `baseml_subtree_summary.tsv`
   - `baseml_tree_manifest.tsv`
   - `overlap_report.tsv`
   - `baseml_subtrees/*.nwk`

### 4.4 拆树校验：`python/validate_phylo_split.py`
该模块会从输入树重新计算期望结果，再与输出文件逐项比对。

它不是简单检查文件存在，而是做“从头复算 + 逐项一致性验证”。

校验包括：
- `core_subtree_summary.tsv` 是否与预期一致
- `baseml_subtree_summary.tsv` 是否与预期一致
- manifest / overlap 是否一致
- 每棵 core tree 的 tip 集和 hash 是否正确
- 每棵 baseml tree 是否 rooted、binary、含唯一 `RSRS`

## 5. 合树模块逻辑
### 5.1 合树主控：`pipe/run_merge_pipeline.sh`
控制流如下：

1. 读取 `Config.json`
2. 检查环境
3. 检查拆树阶段必需文件是否存在：
   - `output/intermediate/rooted.tree`
   - `output/core_subtree_summary.tsv`
   - `output/baseml_subtree_summary.tsv`
   - `output/baseml_tree_manifest.tsv`
4. 预清理 `output/merge/`
5. 若 `analysis_tree_source=simulated`，先调用 `simulate_baseml_results.py`
6. 调用 `merge_baseml_subtrees.py`
7. 调用 `validate_merged_tree.py`

### 5.2 合树公共库：`python/phylo_merge_common.py`
这个模块提供合树阶段所有可复用逻辑。

#### A. TSV 数据结构
- `CoreSummaryRecord`
- `BasemlSummaryRecord`

#### B. 结果树模拟
- `stable_hash_int()`
- `randomize_tree_branch_lengths()`

这里使用对数正态乘性扰动：
- 每条边先保证为正值
- 再乘以 lognormal 倍数

#### C. core clade 识别与匹配
- `find_core_mrca()`
- `build_clade_signature_maps()`

核心思想是：
- 在结果子树中先定位 `core_tip_names` 的 MRCA
- 再用“后代 tip 集的 SHA256 signature”匹配主树中的对应边

这是合树能稳定回填 branch length 的关键。

#### D. scaffold 重标定
- `compute_multiplier()`
- `weighted_geometric_mean()`
- `build_descendant_core_ids()`

逻辑是：

1. 每个 `core subtree` 回填后，都能得到一个 `core root multiplier`
2. 对于没有被某个单一 core 直接覆盖的 scaffold 边：
   - 找到它后代范围内所有 core roots
   - 取这些倍率的按 `core_n_tips` 加权几何均值
   - 用这个尺度去更新 scaffold 边

#### E. 合树合法性
- `validate_tree_against_expected_tip_set()`
- `validate_rooted_binary_tree()`
- `validate_branch_lengths_complete()`

### 5.3 模拟模块：`python/simulate_baseml_results.py`
该模块只做一件事：把现有 `baseml_subtrees/*.nwk` 随机化为“模拟 PAML 输出”。

执行流程：

1. 读取 `baseml_subtree_summary.tsv`
2. 逐棵读取 `baseml subtree`
3. 验证 tip 集和 rooted-binary 结构
4. 对 branch lengths 做随机扰动
5. 输出到 `output/merge/simulated_baseml_subtrees/`
6. 记录 `simulation_manifest.tsv`

### 5.4 合并模块：`python/merge_baseml_subtrees.py`
这个模块是真正的“回填器”。

逻辑分两段：

#### A. core 边精确回填
对每个 `baseml subtree`：

1. 在 master tree 中找到对应的 `core_root_node_id`
2. 在结果树中定位 `core_mrca`
3. 计算 master core 和 result core 的 clade signatures
4. 要求二者 signature 集合完全一致
5. 按 signature 一一匹配后，把所有 core 内部边、tip 边、core root 边的长度回填回 master tree

#### B. scaffold 边启发式重标定
对不在任何单一 core 内部的边：

1. 找到该边后代范围内所有 core roots
2. 取这些 core root 的 branch-length multiplier
3. 按 `core_n_tips` 加权做几何均值
4. 用该倍率乘以 scaffold 原始 branch length

最终输出：
- `merged_tree.nwk`
- `merge_report.tsv`
- `edge_update_report.tsv`

### 5.5 合树校验：`python/validate_merged_tree.py`
该模块验证 merged tree 是否符合设计预期。

检查内容包括：
- tip 数是否与 master tree 一致
- tip 集是否完全一致
- 是否 rooted
- `RSRS` 是否仍为根下 singleton child
- 是否严格 binary
- 拓扑是否与 master tree 完全一致
- 所有非 root 边是否有非负 branch length
- `merge_report.tsv` 是否覆盖全部 core / baseml 子树
- `edge_update_report.tsv` 是否覆盖全部非 root 边

## 6. 关键设计思想
### 6.1 为什么拆成两层
- `core partition`
  - 用于全覆盖、无重叠地给每个原始 tip 找唯一归属
- `baseml subtree`
  - 用于实际运行 `PAML baseml`
  - 允许因锚点和外群而重叠

### 6.2 为什么合树不直接拼接整棵 baseml 子树
因为 `baseml subtree` 含有 anchor 和 `RSRS`，它们不是目标核心片段的唯一所有者。若直接把整棵树拼回主树，会引入重复和结构冲突。

因此当前实现只把每棵子树中对应 `core clade` 的 branch lengths 回填回主树。

### 6.3 为什么 scaffold 边做几何均值缩放
因为 scaffold 边没有唯一的单棵结果树可以直接提供 branch length。当前策略用后代已更新 core roots 的倍率来推断其尺度，是一个兼顾稳定性和可解释性的启发式方案。

## 7. 维护建议
### 如果要支持真实 PAML 输出
只需要：

1. 把 `merge.analysis_tree_source` 设为 `external`
2. 准备目录，命名为 `{baseml_subtree_id}.nwk`
3. 保证每棵结果树：
   - tip 集与对应 `baseml_subtree_summary.tsv` 一致
   - rooted
   - binary
   - `RSRS` 唯一存在

### 如果要支持多分叉树
当前项目不适合直接支持多分叉合树，因为：
- `baseml subtree`
  - 当前被显式要求 rooted binary
- merge
  - 当前也显式要求 master tree 与结果树 binary

若要扩展，必须整体修改：
- rooted/binary 校验函数
- core topology signature 比对逻辑
- `baseml` 输入预处理策略

## 8. 快速索引
### 拆树相关
- 主控：`pipe/run_split_pipeline.sh`
- 主模块：`python/split_phylo_tree.py`
- 公共库：`python/phylo_split_common.py`
- 校验：`python/validate_phylo_split.py`

### 合树相关
- 主控：`pipe/run_merge_pipeline.sh`
- 模拟：`python/simulate_baseml_results.py`
- 合并：`python/merge_baseml_subtrees.py`
- 公共库：`python/phylo_merge_common.py`
- 校验：`python/validate_merged_tree.py`

