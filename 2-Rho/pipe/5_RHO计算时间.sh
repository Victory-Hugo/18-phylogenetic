#!/bin/bash
BASE_DIR='/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/2-Rho/'

/home/luolintao/miniconda3/bin/python3 \
    ${BASE_DIR}/script/5_RHO计算时间.py \
    -i ${BASE_DIR}/output/Rho统计.tsv \
    -o ${BASE_DIR}/output/Rho统计时间.tsv