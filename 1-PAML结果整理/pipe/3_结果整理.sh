#!/bin/bash

PYTHON3="/home/luolintao/miniconda3/bin/python3"

BASE_DIR='/mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/'
RESULT_FILE="${BASE_DIR}/output/Example.txt"


"$PYTHON3" \
    ${BASE_DIR}/pipe/python/3_结果整理.py \
    ${RESULT_FILE} \
    ${BASE_DIR}/output/

"$PYTHON3" \
    ${BASE_DIR}/pipe/python/4_分子钟校准.py \
    -t ${BASE_DIR}/output/ID_Length.tree \
    -c ${BASE_DIR}/output/TIP_Length.csv \
    -m 2.53e-8 \
    -o ${BASE_DIR}/output/ \

