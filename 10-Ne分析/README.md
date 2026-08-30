# Ne 轨迹聚类与可视化
对 BEAST Skygrid 的后验结果做进一步分析。

## 快速开始

```bash
bash pipe/run.sh                        # 用 conf/config.yaml
NE_CONFIG=conf/other.yaml bash pipe/run.sh
NE_PYTHON=/path/to/python3 bash pipe/run.sh   # 指定解释器          
```

## 输入

最小输入是 `input/` 下的 **`*.log` 加同名 `*.xml`**：

- `.log`：BEAST 直出的后验日志，必须含 `skygrid.logPopSize*` 列。一个 `.log` = 一个类别，
  类别名取文件名。
- `.xml`：Skygrid 时间轴的权威来源。`.log` 里只有 logPopSize 的取值、没有任何时间信息，
  必须由 XML 的 `skygrid.cutOff` 与网格数才能还原每一段对应哪个时间区间；两者的维度会
  互相核对，对不上直接报错。

可选的 `input/samples.tsv`（`category / display_name / group / n_tips`）只在需要下面三件事
之一时才必要，缺失不影响主流程：

- `display_name`：图表里的展示名。不给则按"下划线换空格并首字母大写"兜底；
- `group`：类别的上层分组，用于弦图外圈色带与热图行注释。不给则全部归一组；
- `n_tips`：开启 `select.size_adjustment` 时才需要。

## 输出

`output/` 下只有两张核心长表与五类图，每张图都配一份同名 TSV（可直接重绘）：

| 文件 | 内容 |
|---|---|
| `⭐1-Feature-Summary.tsv` | 类别 × 时间窗 × 8 个特征的中位数、95% HPD、方向支持度、覆盖度 |
| `⭐2-Cluster-And-Similarity.tsv` | 簇归属、同簇频率、合并分数与权重、降维坐标、解释率、结点大小 |
| `1-Period-Feature-Heatmap` | 每个时间段一个"类别 × 8 特征"热图子图，全部排进同一张图 |
| `2-Composite-Score-Heatmap` | 类别 × 时间 bin 的水平分数与变化分数两个子图，行按共同可靠时段聚类 |
| `3-Cluster-Scatter` | 每个时间段一个 PCA（默认）或 UMAP 散点子图，按簇着色并加置信椭圆 |
| `4-Chord-Diagram` | 每个时间段一个分层弦图子图，弦的粗细为同簇频率 |
| `5-Pairwise-Matrix-Heatmap` | 每个时间段一对"类别 × 类别"矩阵：左同簇频率、右特征相似度 |

所有图都按 **A4 纵向分页**：PDF 是多页的（内容比 A4 宽时整页等比缩到页宽，需要
`pdfjam` 或 `gs`，两者都没有时页面保持 A4 比例的放大版），PNG 每页一个文件
（`*-p01.png`、`*-p02.png`…，只有一页时不带页码）。每页放几个子图由
`figure.panels_per_row` 与各图自身的子图尺寸共同决定——列数越少，缩放越轻，
A4 上的字号越大。
