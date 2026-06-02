# mtDNA Backbone-Target-Graft Pipeline

本项目采用 **骨架树 + 目标 clade graft** 策略，将大树分拆为多棵 PAML 子树独立运行 baseml，再合并为完整 ML 树，最后做超度量标准化和时间校正。

## 流程概览

```
输入大树 + FASTA + 外群
→ 生成 backbone tree（自动抽样或手动指定）
→ 递归切分 target subtree
→ 构造 PAML 子树（backbone + target + anchors + outgroup）
→ 批量运行 baseml
→ 聚合 backbone 边长 + graft 各 target clade → merged tree
→ 超度量标准化 → 时间校正
```

核心原则：每棵 PAML 子树保留同一套骨架树，target clade 只负责自己的非骨架 tips。

## 配置文件

| 阶段 | 配置 |
|------|------|
| 拆树 | [conf/1-split.yaml](conf/1-split.yaml) |
| PAML | [conf/2-paml.yaml](conf/2-paml.yaml) |
| 合树 | [conf/3-merge.yaml](conf/3-merge.yaml) |
| 超度量 | [conf/4-ultrastandard.yaml](conf/4-ultrastandard.yaml) |
| 时间校正 | [conf/5-time_calib.yaml](conf/5-time_calib.yaml) |

关键默认参数：
- `max_tips = 300`，`backbone_size = 100`
- 骨架抽样策略：`greedy_farthest_patristic`
- PAML 模板：`conf/baseml_mtDNA.ctl`（修改 PAML 参数直接改模板即可）
- 合并模式：`backbone_graft`
- 时间校正默认用分子钟法，速率 `2.67e-8`，序列长度 `16569`

## 输入要求

- **树**：Newick 格式，未定根需提供外群列表
- **FASTA**：header 必须与树 tip 对应
- **外群**：纯文本，每行一个 tip，支持 `#` 注释和空行

## 运行

五个阶段独立运行，也可一键串行：

```bash
bash pipe/0-run_all_pipeline.sh                           # 一键运行 1-5
bash pipe/1-run_split_pipeline.sh                         # 拆树
bash pipe/2-run_paml_pipeline.sh                          # PAML
bash pipe/3-run_merge_pipeline.sh                         # 合树
bash pipe/4-run_ultrastandard_pipeline.sh                 # 超度量
bash pipe/5-run_time_calib_pipeline.sh                    # 时间校正
```

指定自定义配置：

```bash
bash pipe/1-run_split_pipeline.sh --config /path/to/1-split.yaml
```

## 注意事项

- **建议优先手动指定骨架**（`paths.backbone_tip_id_file`），而非完全依赖自动抽样
- 时间校正不改变拓扑，若年代偏差大，优先检查上游超度量树的 root-to-tip 深度和分子钟速率
- 骨架自动抽样策略：先将 ingroup 递归切成 `backbone_size` 个前沿 clade，每 clade 选最深代表 tip，避免碎片化
- target 切分依据是非骨架 tip 数，`target_capacity = max_tips - backbone 数 - anchors - 1(外群)`
