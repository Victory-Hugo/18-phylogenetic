# mtDNA Backbone-Target-Graft Pipeline（BEAST 版）

本项目由 `5-拆树PAML合树` 移植而来：采用 **骨架树 + 目标 clade graft** 策略，把大树拆为多棵子树，
**用 BEAST(v10.5.0) 在固定拓扑上贝叶斯估计枝长/分歧时间**（替换原 PAML baseml 步骤），再合并为完整树，
最后做超度量标准化与时间校正。

## 流程概览

```
输入大树 + FASTA + 外群
→ 生成 backbone tree（自动抽样或手动指定）
→ 递归切分 target subtree
→ 为每棵子树（含 backbone）构造固定拓扑 BEAST 作业
→ 批量运行 beast + treeannotator → MCC 超度量子树
→ 聚合 backbone 边长 + graft 各 target clade → merged tree
→ 超度量标准化 → 时间校正
```

## 配置文件

| 阶段 | 配置 |
|------|------|
| 拆树 | [conf/1-split.yaml](conf/1-split.yaml) |
| BEAST | [conf/2-beast.yaml](conf/2-beast.yaml) |
| 合树 | [conf/3-merge.yaml](conf/3-merge.yaml) |
| 超度量 | [conf/4-ultrastandard.yaml](conf/4-ultrastandard.yaml) |
| 时间校正 | [conf/5-time_calib.yaml](conf/5-time_calib.yaml) |

关键默认参数：
- BEAST 模板：`conf/beast_template.xml`（修改替换模型/时钟/先验/算子直接改模板）
- BEAST 模型：GTR + Gamma(4)，strict clock（速率固定 1.0），Constant-size Coalescent
- MCMC：`chain_length = 5e6`，`burnin_percent = 10`，MCC 节点高度取 `mean`
- 合并模式：`backbone_graft`；时间校正默认节点标定（RSRS 180000 年）

## 环境与依赖

- conda 环境：`BigLin`
- BEAST v10.5.0：`/home/luolintao/opt/BEASTv10.5.0/bin/`（`beast`、`treeannotator`，BigLin PATH 可用）
- `gotree`（NEXUS→Newick、拆树辅助）；Python：Biopython 等

## 输入要求

- **树**：Newick；本仓库用已定根的 `data/mtDNA.reroot.tree`
- **FASTA**：`data/mtDNA.fasta`，header 与树 tip 一一对应（`seq_id_strategy: exact`）
- **外群**：`META/outgroup_tips.tsv`（本仓库为 `RSRS`），每行一个 tip，支持 `#` 注释与空行

## 运行

```bash
bash pipe/0-run_all_pipeline.sh                           # 一键运行 1-5
bash pipe/1-run_split_pipeline.sh                         # 拆树
bash pipe/2-prepare_beast_xml.sh                          # 生成 BEAST XML
bash pipe/2-run_beast_xml.sh                              # 运行 BEAST XML 并解析
bash pipe/3-run_merge_pipeline.sh                         # 合树
bash pipe/4-run_ultrastandard_pipeline.sh                 # 超度量
bash pipe/5-run_time_calib_pipeline.sh                    # 时间校正
```

指定自定义配置：`bash pipe/2-prepare_beast_xml.sh --config /path/to/2-beast.yaml`，
然后运行 `bash pipe/2-run_beast_xml.sh --config /path/to/2-beast.yaml`

### fake 快速测试模式

BEAST MCMC 较慢，联调下游时可把 `conf/2-beast.yaml` 的 `beast.execution_mode` 改为 `fake`，
跳过真实 MCMC、直接由子树拓扑随机化生成超度量树，秒级跑通 stage1→3→4→5 全链路。
