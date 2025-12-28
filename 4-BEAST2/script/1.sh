#!/usr/bin/bash

BEAST2_EXE="/mnt/f/Onedrive/文档（科研）/脚本/Download/18-phylogenetic/4-BEAST2/bin/beast"
WORK_DIR="/mnt/f/Onedrive/文档（科研）/脚本/Download/18-phylogenetic/4-BEAST2/input/"
XML_FILE="/mnt/f/Onedrive/文档（科研）/脚本/Download/18-phylogenetic/4-BEAST2/input/testBSP.xml"
cd ${WORK_DIR}

"${BEAST2_EXE}" \
    -beagle_GPU \
    -threads 16 \
    ${XML_FILE}