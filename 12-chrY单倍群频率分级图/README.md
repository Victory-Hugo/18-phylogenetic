# chrY 单倍群频率分级图

按 Y 染色体单倍群的系统发育层级统计各群体频率，并生成分级堆叠柱状图；同时可通过重复分层交叉验证选择每个群体的自适应最佳分辨率。

## 输入与配置

示例输入为 `input/example.tsv`（制表符分隔），至少包含 `ID`、`Haplogroup`、`Class` 和 `Label` 四列。请在 `conf/1-haplogroup_level_frequency.yaml` 与 `conf/2-optimal_depth.yaml` 中设置输入路径、列名、阈值及 Python/R 路径；颜色配置位于 `conf/color.tsv` 和 `conf/label_color.tsv`；`conf/afm/` 存放 Arial 字体度量，供 PDF 输出注册可编辑的 Arial 文本。

## 运行

```bash
bash pipe/1-haplogroup_level_frequency.sh  # 各固定层级的频率图
bash pipe/2-optimal_depth.sh                # 自适应最佳层级及频率图
```

结果表写入 `output/result/`，图和绘图用频率表写入 `output/figure/`。带 `⭐` 前缀的 TSV 是主要输出文件。

