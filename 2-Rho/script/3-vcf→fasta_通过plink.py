#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from Bio import SeqIO

def create_mapping(mapping_path):
    """
    读取 TSV 文件，构建原始 ID 到修正 ID 的映射字典。
    新 ID = 第一列（去除空白） + "_" + 第二列（去除空白）
    
    参数:
        mapping_path (str): TSV 映射文件路径，要求文件有表头且至少包含两列。
    
    返回:
        dict: 键为原始 ID，值为拼接后的新 ID。
    """
    try:
        # 保证所有数据以字符串方式读取
        df = pd.read_csv(mapping_path, sep='\t', header=0, dtype=str)
    except Exception as e:
        raise RuntimeError(f"读取映射文件 {mapping_path} 失败: {e}")

    if df.shape[1] < 2:
        raise ValueError("映射文件的列数不足2，无法构建映射关系。")

    # 创建新ID：第一列与第二列拼接（去除空白），中间用下划线连接
    df['new_id'] = df.iloc[:, 0].str.strip() + "_" + df.iloc[:, 1].str.strip()
    # 将第一列作为原始ID构建映射字典
    mapping = df.set_index(df.columns[0])['new_id'].to_dict()
    return mapping

def process_fasta_records(input_fasta, mapping):
    """
    生成器函数：遍历 FASTA 文件中的记录，如果记录的 ID 在映射字典中，则更新为新 ID。
    
    参数:
        input_fasta (str): 输入 FASTA 文件路径。
        mapping (dict): 原始 ID 到新 ID 的映射字典。
    
    生成:
        更新后的 SeqRecord 对象。
    """
    with open(input_fasta, 'r', encoding='utf-8') as fin:
        for record in SeqIO.parse(fin, 'fasta'):
            if record.id in mapping:
                new_id = mapping[record.id]
                record.id = new_id
                record.name = new_id         # 保证 name 与 id 一致
                record.description = ""    # 更新描述信息
            yield record

def update_fasta_ids(input_fasta, output_fasta, mapping):
    """
    根据映射关系更新 FASTA 文件中的记录 ID，并写入到新的 FASTA 文件中。
    
    参数:
        input_fasta (str): 输入 FASTA 文件路径。
        output_fasta (str): 输出 FASTA 文件路径。
        mapping (dict): 原始 ID 到新 ID 的映射字典。
    """
    with open(output_fasta, 'w', encoding='utf-8') as fout:
        SeqIO.write(process_fasta_records(input_fasta, mapping), fout, 'fasta')

def main():
    parser = argparse.ArgumentParser(
        description="使用映射文件更新 FASTA 文件中的记录 ID"
    )
    parser.add_argument(
        "--mapping",
        required=True,
        help="TSV 映射文件路径，必须包含至少两列，用于构建新ID（新ID = 第一列 + '_' + 第二列）"
    )
    parser.add_argument(
        "--input_fasta",
        required=True,
        help="输入 FASTA 文件路径"
    )
    parser.add_argument(
        "--output_fasta",
        required=True,
        help="输出 FASTA 文件路径"
    )
    
    args = parser.parse_args()

    # 生成映射字典，不产生中间文件
    mapping = create_mapping(args.mapping)
    
    # 更新 FASTA 文件中的记录 ID
    update_fasta_ids(args.input_fasta, args.output_fasta, mapping)
    
    print(f"FASTA 文件已更新，输出路径: {args.output_fasta}")

if __name__ == "__main__":
    main()
