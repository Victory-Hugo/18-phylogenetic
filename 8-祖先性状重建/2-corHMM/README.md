# corHMM祖先性状重建示例

本项目使用一个207末端的D1a1a最大似然子树，演示离散性状的ER、SYM、ARD模型比较、随机映射以及逐事件树枝定位。

## 输入

- `input/1-example.tree`：有根Newick树。
- `input/2-traits.tsv`：仅含`species`和`state`两列，ID与树末端一一对应。

示例性状包括`YellowRiver`、`GansuQinghai`、`SouthWest`和`Plateau`四个状态。

## 运行

所有可调参数均位于`conf/1-corHMM.yaml`。默认按AIC选择模型，根先验为`maddfitz`，使用12个随机起点并执行200次随机映射。选择指标可改为`BIC`，根先验可改为`yang`或`flat`。模型拟合默认使用1核，以保证同一随机种子可以完全复现；可调高`models.n_cores`加速，但corHMM 2.8的并行随机起点可能产生不完全相同的最优值。

```bash
bash pipe/1-corHMM.sh
```

R脚本通过`/usr/bin/Rscript --vanilla`运行，需要`ape`、`corHMM`和`knitr`。
`flat`在模型拟合中直接传给`corHMM()`；随机映射阶段使用数学等价的标量等权重，以规避corHMM 2.8中`makeSimmap(root.p="flat")`的已知实现错误。

## 输出

### `output/result/1-model-compare`

- `1-model-comparison.tsv`：模型参数数、logLik、AIC、BIC、差值和排名。
- `1-best-model.tsv`：最终选中的模型。
- `1-best-model-rate-matrix.tsv`：最优模型的连续时间马尔可夫速率矩阵长表。

### `output/result/2-simmap`

- `2-simmap-transition-summary.tsv`：方向性转移次数的均值、标准差、95%区间和非零比例。
- `2-simmap-edge-support.tsv`：每条父子枝上每种方向性转移的模拟支持率。
- `2-simmap-maps.rds`：完整随机映射对象，可用于复核或后续分析。

### `output/result/3-trans`

- `3-transition-events.tsv`：每次模拟中每个转移事件的父子节点、树枝、方向及枝内位置。
- `3-tree-edges.tsv`：分析树的edge ID、父子节点、枝长和代表性后代末端。
- `3-analysis.tree`：实际用于拟合和随机映射的完全二分分析树。

逐事件位置以树的枝长为单位，并非日历年代或地理坐标。`relative_position_from_parent=0`表示靠近父节点，`1`表示靠近子节点。输入树的多分叉会以确定性方式解析；由此产生的零长枝会按配置设置下限，并在结果中以`edge_was_floored`标记。
