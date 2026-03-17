library(ggtree)

tree <- read.tree("/mnt/f/onedrive/文档（科研）/脚本/Download/18-phylogenetic/5-拆树PAML合树/script/tmp.tree")

tree  |> ggtree()
