# -*- coding: utf-8 -*-
import sys
import re
import csv
from pathlib import Path
import pandas as pd

# ===================== 基础提取功能 =====================
def extract_after_keyword(input_path: Path, keyword: str, offsets: list[int], output_paths: list[Path]):
    lines = input_path.read_text(encoding='utf-8').splitlines()
    for idx, line in enumerate(lines):
        if keyword in line:
            for offset, out_path in zip(offsets, output_paths):
                if idx + offset < len(lines):
                    out_path.write_text(lines[idx + offset] + '\n', encoding='utf-8')

def extract_multiple_lines_after_keyword(input_path: Path, keyword: str, next_lines: list[int], output_path: Path):
    lines = input_path.read_text(encoding='utf-8').splitlines()
    extracted = []
    for idx, line in enumerate(lines):
        if keyword in line:
            for offset in next_lines:
                if idx + offset < len(lines):
                    extracted.append(lines[idx + offset])
    output_path.write_text('\n'.join(extracted) + '\n', encoding='utf-8')

# ===================== 数据处理功能 =====================
def read_branch_data(parent_child_path: Path, se_path: Path) -> pd.DataFrame:
    lines = [l.strip() for l in parent_child_path.read_text(encoding='utf-8').splitlines() if l.strip()]
    pcs = lines[0].split()
    params = lines[1].split()
    ses = se_path.read_text(encoding='utf-8').split()
    n = min(len(pcs), len(params), len(ses))
    return pd.DataFrame([
        {
            'Branch': tag,
            'Param': param,
            'SE': se,
            'Parent': int(tag.split('..')[0]),
            'Child': int(tag.split('..')[1])
        }
        for tag, param, se in zip(pcs[:n], params[:n], ses[:n])
    ])

def read_tip_mapping(tree_path: Path) -> pd.DataFrame:
    text = tree_path.read_text(encoding='utf-8')
    matches = re.findall(r'(\d+)_([^,\)\s]+)', text)
    df = pd.DataFrame(matches, columns=['Node', 'ID'])
    df['Node'] = df['Node'].astype(int)
    return df

def merge_branch_and_tips(branch_df: pd.DataFrame, tip_df: pd.DataFrame) -> pd.DataFrame:
    return branch_df.merge(
        tip_df.rename(columns={'Node': 'Parent', 'ID': 'Parent_ID'}), on='Parent', how='left'
    ).merge(
        tip_df.rename(columns={'Node': 'Child', 'ID': 'Child_ID'}), on='Child', how='left'
    )

def extract_tip_lengths(tree_path: Path, output_csv: Path):
    text = tree_path.read_text(encoding='utf-8')
    matches = re.findall(r'([A-Za-z0-9\-_]+):\s*([0-9]*\.[0-9]+)', text)
    with output_csv.open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['ID', 'Length'])
        writer.writerows(matches)

def clean_temp_files(files: list[Path]):
    for file in files:
        if file.exists():
            file.unlink()
            print(f"已删除临时文件：{file.name}")

# ===================== 主程序 =====================
def main():
    if len(sys.argv) != 3:
        print("用法: python script.py <输入文件路径> <输出目录路径>")
        sys.exit(1)

    input_path = Path(sys.argv[1]).resolve()
    output_dir = Path(sys.argv[2]).resolve()

    if not input_path.is_file():
        print(f"错误：输入文件不存在：{input_path}")
        sys.exit(1)
    if not output_dir.is_dir():
        print(f"错误：输出目录不存在：{output_dir}")
        sys.exit(1)

    # 定义所有输出文件
    parent_child_file = output_dir / 'PC.txt'
    se_file = output_dir / 'SE.txt'
    tip_id_file = output_dir / 'TIP_ID.tree.txt'

    id_length_file = output_dir / 'ID_Length.tree'
    merged_csv = output_dir / 'Parameters.csv'
    tip_length_csv = output_dir / 'TIP_Length.csv'

    temp_files = [parent_child_file, se_file, tip_id_file]

    # Step 1. 提取中间文件
    extract_multiple_lines_after_keyword(input_path, 'lnL', [1, 2], parent_child_file)
    extract_multiple_lines_after_keyword(input_path, 'SEs for parameters', [1], se_file)
    extract_after_keyword(input_path, 'tree length', [4, 6], [id_length_file, tip_id_file])

    # Step 2. 处理数据生成最终文件
    branch_df = read_branch_data(parent_child_file, se_file)
    tip_df = read_tip_mapping(tip_id_file)
    merged_df = merge_branch_and_tips(branch_df, tip_df)
    merged_df.to_csv(merged_csv, index=False, encoding='utf-8-sig')
    print(f"已生成合并表格：{merged_csv}")

    extract_tip_lengths(id_length_file, tip_length_csv)
    print(f"已生成 tip 长度表：{tip_length_csv}")

    # Step 3. 删除临时文件
    clean_temp_files(temp_files)

if __name__ == '__main__':
    main()
