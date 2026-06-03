#!/usr/bin/env python3
"""骨架(backbone)树校准定年作业准备 —— 生成"校准一致热启动"的单个 BEAST XML。

读取 conf/2-beast.yaml 的 ``beast.calibration`` 配置，对骨架树：
  1. 用 phylotree v17.2 lineage 判定 M/N 单系成员并校验单系性；
  2. 复用 python/beast_common.py 生成 taxa / alignment / 超度量起始树（substitutions）；
  3. 线性预缩放到「年」(root=root.age)，并把 M/N clade 重置到各自校准年龄 → 校准一致**热启动**起始树；
  4. 注入 UCLD 宽松钟 + root/M/N 三点 normalPrior 校准的模板；
  5. 输出 1 个 XML（file_prefix=backbone_only），与正式 stage2 backbone 作业保持一致。

不修改正式 split/merge/stage5 逻辑；仅复用 beast_common 的工具函数。

用法:
    python3 python/prepare_calibrated_backbone_job.py [--config conf/2-beast.yaml]
"""
from __future__ import annotations

import argparse
import io
import json
import sys
from pathlib import Path

import yaml
from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent))

from beast_common import (
    build_alignment_block,
    build_taxa_block,
    build_ultrametric_starting_newick,
    map_tip_to_source_id,
    read_fasta_record_map,
)

PLACEHOLDER_UCLD_STDEV_START = "0.1"


# --------------------------------------------------------------------------- #
# M/N 归属判定（phylotree lineage）与单系校验
# --------------------------------------------------------------------------- #
def classify_M_N(tip_names, id_hap_tsv: Path, phylotree_json: Path):
    id2hap = {}
    for line in Path(id_hap_tsv).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        id2hap[parts[0]] = parts[1]
    J = json.loads(Path(phylotree_json).read_text(encoding="utf-8"))["haplogroups"]

    def lineage(hap):
        if not hap:
            return None
        if hap in J:
            return J[hap]["lineage"]
        for cut in range(len(hap), 0, -1):  # 逐级回退尾部修饰
            if hap[:cut] in J:
                return J[hap[:cut]]["lineage"]
        return None

    M, N, miss = [], [], []
    for tp in tip_names:
        if tp == "RSRS":
            continue
        lin = lineage(id2hap.get(tp))
        if lin is None:
            miss.append(tp)
        elif "M" in lin:
            M.append(tp)
        elif "N" in lin:
            N.append(tp)
    return M, N, miss


def check_monophyly(tree, members, label):
    anc = tree.common_ancestor(members)
    desc = {t.name for t in anc.get_terminals()}
    extra = desc - set(members) - {"RSRS"}
    ok = extra == set()
    print(f"[mono] {label}: 成员={len(members)} MRCA后代={len(desc)} 混入={len(extra)} 单系={ok}")
    if extra:
        print(f"        混入样例: {sorted(extra)[:8]}")
    return ok


# --------------------------------------------------------------------------- #
# 起始树：超度量(substitutions) → 预缩放到年 → M/N clade 重置到校准年龄（热启动）
# --------------------------------------------------------------------------- #
def _node_heights(tree):
    dep = {}

    def rec(c, d):
        dep[c] = d
        for ch in c.clades:
            rec(ch, d + (ch.branch_length or 0.0))

    rec(tree.root, 0.0)
    rooth = max(dep[t] for t in tree.get_terminals())
    return dep, rooth


def _rescale_clade(tree, members, target, label):
    anc = tree.common_ancestor(members)
    dep, rooth = _node_heights(tree)
    cur = rooth - dep[anc]
    factor = target / cur
    for c in anc.find_clades():
        if c is anc:
            continue
        if c.branch_length:
            c.branch_length *= factor
    anc.branch_length = (anc.branch_length or 0.0) + (cur - target)  # 茎枝补偿 → anc 高度=target
    print(f"[warm] {label}: 起始高度 {cur:.0f} -> {target:.0f} (factor={factor:.3f})")


def build_warmstart_newick(subst_newick, M, N, root_age, M_age, N_age, warm_start: bool):
    """返回 (warmstart_newick_years, ucld_mean_start, root_subst)。"""
    tree = Phylo.read(io.StringIO(subst_newick), "newick")
    _, root_subst = _node_heights(tree)
    scale = root_age / root_subst
    for c in tree.find_clades():
        if c.branch_length is not None:
            c.branch_length *= scale  # → 年尺度，root=root_age
    if warm_start:
        _rescale_clade(tree, M, M_age, "M")
        _rescale_clade(tree, N, N_age, "N")
    else:
        print("[warm] warm_start=false：使用 IQTREE 压缩起点（不推荐，易被钉死）")
    buf = io.StringIO()
    Phylo.write(tree, buf, "newick", format_branch_length="%.6f")
    newick = buf.getvalue().strip().replace("\n", "")
    return newick, root_subst / root_age, root_subst


# --------------------------------------------------------------------------- #
# 校准统计量块
# --------------------------------------------------------------------------- #
def build_calibration_block(M, N):
    def taxa(name, tips):
        rows = "\n".join(f'\t\t<taxon idref="{t}"/>' for t in tips)
        return f'\t<taxa id="{name}">\n{rows}\n\t</taxa>'

    def tmrca(stat_id, taxa_id):
        return (f'\t<tmrcaStatistic id="{stat_id}" absolute="true">\n'
                f'\t\t<mrca><taxa idref="{taxa_id}"/></mrca>\n'
                f'\t\t<treeModel idref="treeModel"/>\n'
                f'\t</tmrcaStatistic>')

    return "\n".join([taxa("M_clade", M), taxa("N_clade", N),
                      tmrca("tmrca_M", "M_clade"), tmrca("tmrca_N", "N_clade")])


# --------------------------------------------------------------------------- #
def run(config_path: Path):
    root_dir = config_path.resolve().parents[1]
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    calib = cfg["beast"]["calibration"]
    if not calib.get("enabled", False):
        print("[skip] beast.calibration.enabled=false")
        return 0

    def rp(p):  # 相对项目根解析
        return (root_dir / p).resolve()

    backbone_tree = rp(calib["backbone_tree"])
    template = rp(calib["template"])
    out_dir = rp(calib.get("output_dir", cfg["beast"].get("backbone_job_dir", "output/beast/backbone_job")))
    out_dir.mkdir(parents=True, exist_ok=True)

    nodes = calib["nodes"]
    root_age, root_sd = float(nodes["root"]["age"]), float(nodes["root"]["stdev"])
    M_age, M_sd = float(nodes["M"]["age"]), float(nodes["M"]["stdev"])
    N_age, N_sd = float(nodes["N"]["age"]), float(nodes["N"]["stdev"])
    n_chains = int(calib.get("n_independent_chains", 1))
    if n_chains != 1:
        print(f"[warn] n_independent_chains={n_chains} 已废弃；按每棵树 1 个 XML 执行。")
        n_chains = 1
    chain_length = int(calib.get("chain_length", cfg["beast"]["chain_length"]))
    log_every = int(calib.get("log_every", cfg["beast"]["log_every"]))
    popsize = float(calib.get("popsize_start", 100000.0))

    # 1. 读骨架树、判定 M/N
    tree = Phylo.read(str(backbone_tree), "newick")
    tip_names = [t.name for t in tree.get_terminals()]
    print(f"[info] 骨架 tips={len(tip_names)} 含RSRS={'RSRS' in tip_names}")
    M, N, miss = classify_M_N(tip_names, rp(calib["id_hap_tsv"]), rp(calib["phylotree_json"]))
    print(f"[info] M={len(M)} N={len(N)} 未匹配lineage={len(miss)}")
    if miss:
        print(f"[warn] 未匹配 lineage: {miss[:8]}")
    okM = check_monophyly(tree, M, "M")
    okN = check_monophyly(tree, N, "N")
    if not (okM and okN):
        print("[warn] M/N 非严格单系：tmrcaStatistic 仍用 MRCA，但校准会纳入混入 tip，结果需谨慎。")

    # 2. taxa / alignment / 超度量起始树（substitutions）
    print("[info] 读取 FASTA（约 163MB，请稍候）…")
    fasta_map = read_fasta_record_map(rp(cfg["paths"]["input_fasta"]))
    tip2src = map_tip_to_source_id(tip_names, fasta_map.keys(), cfg["beast"].get("seq_id_strategy", "exact"))
    taxa_block = build_taxa_block(tip_names)
    align_block = build_alignment_block(tip_names, fasta_map, tip2src)
    subst_newick = build_ultrametric_starting_newick(
        backbone_tree, len(tip_names), min_branch_length=1e-8, ultrametric_tolerance=1e-8)

    # 3. 热启动起始树（年）
    warm_newick, ucld_mean, root_subst = build_warmstart_newick(
        subst_newick, M, N, root_age, M_age, N_age, calib.get("warm_start", True))
    print(f"[info] root_subst={root_subst:.6g} ucld.mean起始={ucld_mean:.6g}")

    calib_block = build_calibration_block(M, N)
    tpl_text = template.read_text(encoding="utf-8")

    # 4. 渲染单个 XML
    written = []
    for i in range(1, n_chains + 1):
        prefix = "backbone_only"
        xml = (tpl_text
               .replace("{{TAXA}}", taxa_block)
               .replace("{{ALIGNMENT}}", align_block)
               .replace("{{STARTING_TREE}}", warm_newick)
               .replace("{{CALIBRATION_BLOCK}}", calib_block)
               .replace("{{UCLD_MEAN}}", f"{ucld_mean:.8g}")
               .replace("{{POPSIZE}}", f"{popsize:.6g}")
               .replace("{{ROOT_AGE}}", f"{root_age:.6g}")
               .replace("{{ROOT_SD}}", f"{root_sd:.6g}")
               .replace("{{M_AGE}}", f"{M_age:.6g}")
               .replace("{{M_SD}}", f"{M_sd:.6g}")
               .replace("{{N_AGE}}", f"{N_age:.6g}")
               .replace("{{N_SD}}", f"{N_sd:.6g}")
               .replace("{{FILE_PREFIX}}", prefix)
               .replace("{{CHAIN_LENGTH}}", str(chain_length))
               .replace("{{LOG_EVERY}}", str(log_every)))
        leftover = [tok for tok in ("{{", "}}") if tok in xml]
        if leftover:
            raise SystemExit(f"[error] XML 仍有未替换占位符: {prefix}")
        dest = out_dir / f"{prefix}.xml"
        dest.write_text(xml, encoding="utf-8")
        written.append(dest)

    # 同时落盘热启动起始树与 M/N tip 清单，便于核查
    (out_dir / "warmstart_starting_tree.years.nwk").write_text(warm_newick + "\n", encoding="utf-8")
    (out_dir / "M_tips.txt").write_text("\n".join(M), encoding="utf-8")
    (out_dir / "N_tips.txt").write_text("\n".join(N), encoding="utf-8")

    print(f"\n[ok] 生成 {len(written)} 个 backbone XML -> {out_dir}")
    for d in written:
        print(f"      {d.name}  ({d.stat().st_size/1e6:.1f} MB)")
    print(f"[run] 远端各链建议: beast -working -overwrite -beagle_SSE -seed <不同种子> {written[0].name}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description="Prepare calibrated backbone BEAST jobs (warm-start, multi-chain).")
    parser.add_argument("--config", default=str(Path(__file__).resolve().parents[1] / "conf/2-beast.yaml"))
    args = parser.parse_args(argv)
    return run(Path(args.config))


if __name__ == "__main__":
    raise SystemExit(main())
