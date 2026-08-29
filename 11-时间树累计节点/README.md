# 时间树累计节点山脊图

把一棵时间树的内部节点按样本分组摊开，估计各组**合并事件随时间的密度**，画成山脊图。
换一棵树 + 一张 metadata TSV，改 `conf/1-conf.yaml` 即可复用。

## 运行

```bash
bash pipe/1-pipe.sh                 # 用默认配置
bash pipe/1-pipe.sh conf/其他.yaml   # 指定配置
```

## 输入

放进 `input/`：

| 文件 | 要求 |
| --- | --- |
| 时间树 | **单棵** newick，枝长单位为年，超度量（ultrametric）。含多棵树会直接报错 |
| metadata | 制表符分隔，一行一个样本，含样本编号列与用于分组的类别列 |

树的 tip 名必须就是样本编号，与 metadata 的 `cohort.id_column` 列取值一致。
实际分析的样本是**两者的交集**：树上没有的 metadata 行、metadata 里没有的 tip 都会被跳过，
数量写在日志里。
