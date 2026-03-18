#!/usr/bin/env python3
import os
import argparse
from Bio import Phylo
import pandas as pd

def convert_tree(tree_in, tree_out, mu):
    tree = Phylo.read(tree_in, "newick")
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length = clade.branch_length / mu
    Phylo.write(tree, tree_out, "newick")
    print(f"[√] 树文件已转换并保存：{tree_out}")

def convert_csv(csv_in, csv_out, mu):
    df = pd.read_csv(csv_in)
    df['Time_years'] = df['Length'] / mu
    df.to_csv(csv_out, index=False)
    print(f"[√] CSV 文件已转换并保存：{csv_out}")

def main():
    parser = argparse.ArgumentParser(
        description="将 Newick 树和 TIP_Length.csv 的分支长 (subs/site) 按突变率换算为历时 (years)"
    )
    parser.add_argument("-t", "--tree",
                        required=True,
                        help="输入 Newick 树文件路径（含分支长），例如 ID_Length.tree")
    parser.add_argument("-c", "--csv",
                        required=True,
                        help="输入 CSV 文件路径（含 ID,Length 列），例如 TIP_Length.csv")
    parser.add_argument("-m", "--mu",
                        type=float,
                        required=True,
                        help="突变率（subs per site per year），如 2.53e-8")
    parser.add_argument("-o", "--outdir",
                        required=True,
                        help="输出目录（若不存在会创建），结果文件保存在此目录下")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    tree_out = os.path.join(args.outdir,
                            os.path.basename(args.tree).replace("Length.tree", "Time.tree"))
    csv_out  = os.path.join(args.outdir,
                            os.path.basename(args.csv).replace("Length.csv", "Time.csv"))

    convert_tree(args.tree, tree_out, args.mu)
    convert_csv(args.csv, csv_out, args.mu)

if __name__ == "__main__":
    main()
