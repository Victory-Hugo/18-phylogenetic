"""
s04_extract_genotype_matrix.py
从VCF提取SNP基因型矩阵，保存为numpy压缩格式。

实现策略：使用 bcftools query 批量输出所有样本的GT，再用numpy解析。
与pysam逐样本循环相比，速度提升约10-20倍。

编码方式：
    0  = REF
    1  = ALT1
    2  = ALT2
    -1 = 缺失（./.）或其他无效基因型

用法（CLI）：
    python s04_extract_genotype_matrix.py \
        --vcf data/merged_clean.vcf.gz \
        --site-masks output/stage_1/intermediate/site_masks.json \
        --output-matrix output/stage_1/intermediate/geno_matrix.npz \
        --output-positions output/stage_1/intermediate/geno_positions.npy \
        --output-samples output/stage_1/intermediate/geno_samples.txt \
        --bcftools /path/to/bcftools

用法（import）：
    from s04_extract_genotype_matrix import run
    run(...)
"""

import argparse
import json
import logging
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pysam

log = logging.getLogger(__name__)

MISSING = np.int8(-1)


def _get_samples_from_vcf(vcf_path: str, bcftools: str) -> list[str]:
    """用bcftools获取VCF中的样本顺序列表。"""
    result = subprocess.run(
        [bcftools, "query", "-l", vcf_path],
        capture_output=True, text=True, check=True,
    )
    return [s.strip() for s in result.stdout.strip().split("\n") if s.strip()]


def _write_sites_file(snp_positions: list[int], chrom: str, tmp_path: str):
    """将SNP位点列表写为bcftools -T 格式（CHROM\tPOS，1-based）。"""
    with open(tmp_path, "w") as f:
        for pos in sorted(snp_positions):
            f.write(f"{chrom}\t{pos}\n")


def _detect_chrom(vcf_path: str, bcftools: str) -> str:
    """从VCF中检测染色体名称（mtDNA通常为 chrM 或 MT）。"""
    result = subprocess.run(
        [bcftools, "view", "-H", "--no-version", "-r", "chrM,MT,NC_012920.1", vcf_path],
        capture_output=True, text=True,
    )
    for line in result.stdout.split("\n"):
        if line.strip() and not line.startswith("#"):
            return line.split("\t")[0]
    # 备用：取第一条记录的CHROM
    result2 = subprocess.run(
        [bcftools, "view", "-H", "--no-version", vcf_path],
        capture_output=True, text=True,
    )
    for line in result2.stdout.split("\n"):
        if line.strip():
            return line.split("\t")[0]
    return "chrM"


def _extract_matrix_bcftools(
    vcf_path: str,
    snp_positions: list[int],
    bcftools: str,
) -> tuple:
    """
    用bcftools query批量提取所有SNP位点的基因型，
    返回 (matrix, ordered_positions, samples)。

    matrix: numpy int8数组，shape=(n_samples, n_snp_sites)
    ordered_positions: POS列表（与matrix列顺序对应，来自bcftools实际输出顺序）
    samples: list[str]，样本ID列表，顺序与matrix行对应
    """
    samples = _get_samples_from_vcf(vcf_path, bcftools)
    n_samples = len(samples)
    log.info(f"  VCF样本数: {n_samples}")

    chrom = _detect_chrom(vcf_path, bcftools)
    log.info(f"  mtDNA染色体名: {chrom}")

    # 写临时位点文件
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".sites.txt", delete=False
    ) as tmp:
        sites_file = tmp.name
        for pos in sorted(snp_positions):
            tmp.write(f"{chrom}\t{pos}\n")

    # 同时输出POS和GT：格式 %POS\t[%GT\t]\n
    # %POS 是该位点的1-based坐标；[%GT] 对每个样本输出GT
    query_fmt = "%POS\\t[%GT\\t]\\n"
    cmd = [
        bcftools, "query",
        "-f", query_fmt,
        "-T", sites_file,
        "-i", 'TYPE="snp"',   # 只保留SNP记录，排除indel
        vcf_path,
    ]
    log.info("运行bcftools query提取基因型矩阵...")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # 逐行解析：每行为一个SNP位点
    positions_out = []
    geno_cols = []    # 每列对应一个位点，shape=(n_samples,)

    n_processed = 0
    for raw_line in proc.stdout:
        line = raw_line.decode("ascii", errors="replace").rstrip("\n")
        if not line:
            continue

        parts = line.split("\t")
        pos = int(parts[0])
        gts = parts[1:]   # 每个样本的GT字符串（如 "0/0", "1/1", ".", "./."）

        # 解析GT：取第一个等位基因字符
        # GT格式：0/0 → '0'→0；1/1 → '1'→1；./. 或 . → '.'→-1
        col = np.empty(n_samples, dtype=np.int8)
        for i, gt_str in enumerate(gts[:n_samples]):
            gt_str = gt_str.strip()
            if not gt_str or gt_str[0] == ".":
                col[i] = MISSING
            else:
                try:
                    allele = int(gt_str[0])
                    col[i] = np.int8(allele) if allele <= 2 else MISSING
                except (ValueError, IndexError):
                    col[i] = MISSING

        positions_out.append(pos)
        geno_cols.append(col)
        n_processed += 1

        if n_processed % 2000 == 0:
            log.info(f"  已处理 {n_processed} 个SNP位点...")

    stdout_remainder, stderr = proc.communicate()
    if proc.returncode != 0:
        log.error(f"bcftools query 报错: {stderr.decode()}")
        raise RuntimeError(f"bcftools query 失败，返回码: {proc.returncode}")

    log.info(f"  bcftools共输出 {n_processed} 个SNP位点")

    if not geno_cols:
        raise ValueError("未提取到任何SNP位点，请检查site_masks和VCF")

    # 矩阵: (n_samples, n_sites)
    matrix = np.stack(geno_cols, axis=1)
    positions = np.array(positions_out, dtype=np.int32)

    Path(sites_file).unlink(missing_ok=True)
    log.info(f"  矩阵 shape: {matrix.shape}，内存: {matrix.nbytes/1e9:.2f} GB")
    return matrix, positions, samples


def run(
    vcf: str,
    site_masks: str,
    output_matrix: str,
    output_positions: str,
    output_samples: str,
    bcftools: str = "bcftools",
) -> int:
    """
    提取VCF基因型矩阵并保存。

    参数：
        vcf:              VCF文件路径
        site_masks:       s01输出的site_masks.json路径
        output_matrix:    输出numpy矩阵路径（.npz）
        output_positions: 输出POS数组路径（.npy）
        output_samples:   输出样本ID列表路径（.txt）
        bcftools:         bcftools可执行路径

    返回：
        0 表示成功
    """
    for path in [output_matrix, output_positions, output_samples]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    log.info(f"读取位点掩码: {site_masks}")
    with open(site_masks) as f:
        masks = json.load(f)
    snp_positions = masks["all_snp_positions"]
    log.info(f"  目标SNP位点数: {len(snp_positions)}")

    log.info(f"从VCF提取基因型矩阵: {vcf}")
    matrix, positions, samples = _extract_matrix_bcftools(vcf, snp_positions, bcftools)

    np.savez_compressed(output_matrix, matrix=matrix)
    log.info(f"基因型矩阵 -> {output_matrix}")

    np.save(output_positions, positions)
    log.info(f"POS数组 -> {output_positions}")

    Path(output_samples).write_text("\n".join(samples) + "\n")
    log.info(f"样本列表 -> {output_samples}")

    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="从VCF提取基因型矩阵（s04）")
    p.add_argument("--vcf",              required=True, help="VCF文件路径")
    p.add_argument("--site-masks",       required=True, help="site_masks.json路径")
    p.add_argument("--output-matrix",    required=True, help="输出矩阵路径（.npz）")
    p.add_argument("--output-positions", required=True, help="输出POS数组路径（.npy）")
    p.add_argument("--output-samples",   required=True, help="输出样本列表路径（.txt）")
    p.add_argument("--bcftools",         default="bcftools", help="bcftools可执行路径")
    return p


def main(argv=None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(
        vcf=args.vcf,
        site_masks=args.site_masks,
        output_matrix=args.output_matrix,
        output_positions=args.output_positions,
        output_samples=args.output_samples,
        bcftools=args.bcftools,
    )


if __name__ == "__main__":
    raise SystemExit(main())
