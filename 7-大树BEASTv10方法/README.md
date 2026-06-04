# mtDNA 大树 BEAST v10 分析流程

本目录用于把人类 mtDNA 全长序列和 IQ-TREE ML 树转换为 BEAST v10 可运行的 Thorney XML。流程主要包括：

1. 生成校准节点成员表。
2. 用 IQ-TREE 3 内置 LSD2 生成时间尺度 starting tree。
3. 生成 BEAST v10 XML，并用短链实跑校验 XML 是否可解析。

## 注意

若你需要使用**chrY**的序列，请**自行准备**一个tsv文件放到`7-大树BEASTv10方法/output/0-calibration_membership/1-table/calibration_membership.tsv`并跳过pipe的第一步：`7-大树BEASTv10方法/conf/0-calibration_membership.yaml`。

`calibration_membership.tsv`的格式是：

```sh
node	sample_id
M	AF346972
M	AF382013
M	AP008264
```

## 输入文件

输入文件放在 `input/`：

| 文件 | 说明 |
|---|---|
| `input/*.fasta` | 序列 |
| `input/*.IQTREE.tree` | IQ-TREE 3 构建并重根后的 ML 树 |


## 校准节点

校准节点写在 `conf/*.yaml` 的 `calibration.nodes` 中。目前使用：

| 节点 | 年代 | 说明 |
|---|---:|---|
| `root` | 180000 | 全树根年龄 |
| `M` | 49000 | 单倍群 M 的 TMRCA |
| `N` | 59000 | 单倍群 N 的 TMRCA |

## 运行流程

在项目根目录运行：

```bash
bash pipe/0-calibration_membership.sh
bash pipe/1-lsd2.sh
bash pipe/2-beast_xml_skygrid.sh
```

如果需要生成 Bayesian Skyline 版本 XML，运行：

```bash
bash pipe/2-beast_xml_bayesian_skyline.sh
```

## 主要输出

### Step 0：校准成员表

```text
output/0-calibration_membership/1-table/calibration_membership.tsv
output/0-calibration_membership/3-report/calibration_membership_report.md
```

### Step 1：LSD2 时间树

```text
output/1-lsd2/1-table/startingTree.lsd2.nwk
output/1-lsd2/1-table/lsd.timetree.nex
output/1-lsd2/1-table/input_validation.tsv
output/1-lsd2/3-report/lsd2_report.md
```

### Step 2：BEAST XML

Skygrid 版本：

```text
output/2-beast_skygrid_xml/1-table/mtDNA.beast.v10.thorney_skygrid.xml
output/2-beast_skygrid_xml/1-table/mtDNA.beast.v10.thorney_skygrid.trees
output/2-beast_skygrid_xml/3-report/build_beast_xml_report.md
```

Bayesian Skyline 版本：

```text
output/2-beast_xml_bayesian_skyline/1-table/mtDNA.beast.v10.thorney.bayesian_skyline.xml
output/2-beast_xml_bayesian_skyline/1-table/mtDNA.beast.v10.thorney.bayesian_skyline.trees
output/2-beast_xml_bayesian_skyline/3-report/build_beast_xml_report.md
```

## 运行 BEAST

生成 Skygrid XML 后，可运行：

```bash
bash script/0-run_beast.sh
```

该脚本会在 `output/2-beast_skygrid_xml/1-table/` 中运行 BEAST，并自动检测 `mtDNA.state_*` 断点文件进行续跑。

如需对同一个 XML 同时启动多个随机种子任务，可使用：

```bash
bash script/0-run_beast.sh \
  --xml output/2-beast_skygrid_xml/1-table/mtDNA.beast.v10.thorney_skygrid.xml \
  --runs 3 \
  --threads 8
```

批量任务会写入 `output/2-beast_skygrid_xml/1-table/jobs/<时间戳>/seed_<随机种子>/`。每个 `seed_*` 目录会复制一份 XML，并独立生成 `mtDNA.beast.stdout.log`、BEAST `.log/.trees` 与 `mtDNA.state_*` 断点文件；复用同一个 `--timestamp` 可让对应 job 目录继续从自己的最新断点续跑。

查看参数：

```bash
bash script/0-run_beast.sh --help
```

续跑后如需合并日志和树文件：

```bash
bash script/1-merge_result.sh
```

若旧单任务目录中没有 `.log/.trees`，脚本会自动选择 `output/2-beast_skygrid_xml/1-table/jobs/` 下最新的时间戳目录进行合并。

批量任务完成后，可合并指定时间戳目录下所有 `seed_*` 的结果：

```bash
bash script/1-merge_result.sh \
  --batch-dir output/2-beast_skygrid_xml/1-table/jobs/<时间戳>
```

批量模式默认把 `merge.log` 和 `merge.tree` 写回该 `<时间戳>` 目录；如需指定输出目录，可增加 `--outdir <目录>`。

合并结果：

```text
output/2-beast_skygrid_xml/merge.log
output/2-beast_skygrid_xml/merge.tree
```
