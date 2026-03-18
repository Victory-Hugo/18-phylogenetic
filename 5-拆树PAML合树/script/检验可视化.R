library(ggtree)

tree <- read.tree("/mnt/l/22-WHALE/8-优化PAML大于ML进行校正方案/output/ultrastandard/merged_ultrametric_tree.nwk")

tree  |> ggtree()
