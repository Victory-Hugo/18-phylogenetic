#!/bin/bash

INPUT_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/output/placer.aln.snp-sites.raxml.jplace"
SCRIPT_R="/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/script/2-可视化.R"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/output/"

Rscript \
    $SCRIPT_R \
    $INPUT_FILE \
    "" \
    $OUTPUT_DIR