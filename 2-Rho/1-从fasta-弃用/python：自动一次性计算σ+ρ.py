import os
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq  # 导入Seq类
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path
import math
import pandas as pd
import numpy as np
import shutil

# 定义桌面路径
DESKTOP_PATH = os.path.join(os.path.expanduser("~"), 'Desktop')
HAPLOGROUP_FILE_PATH = os.path.join('F:/OneDrive/文档（科研）/脚本/我的科研脚本/Python/母系专用/线粒体单倍群phylotree(version17)2024年8月19日.txt')
INPUT_TXT_PATH = r'C:/Users/victo/Desktop/新建 Text Document.txt'
FASTA_FILE_PATH = r'C:/Users/victo/Desktop/古代现代对齐.fasta'
OUTPUT_FASTA_PATH = None  # 根据用户输入的单倍群动态生成

# 输出文件路径
OUTPUT_TXT_FILE = Path('C:/Users/victo/Desktop/output_results.txt')
FINAL_OUTPUT_FILE = Path(f'C:/Users/victo/Desktop/Rho统计共祖时间.txt')
RHO_RESULT_FILE = Path('C:/Users/victo/Desktop/ρ结果.txt')

# ------------------- 单倍群提取部分 -------------------
def parse_haplogroup_file(lines):
    print("正在解析Haplogroup文件...")
    haplogroups = []
    for line in lines:
        stripped_line = line.strip()
        if stripped_line:
            level = line.count('\t')
            haplogroups.append((level, stripped_line))
    return haplogroups

def find_downstream_haplogroups(haplogroup_name, haplogroups):
    print("正在查找下级单倍群...")
    downstream_haplogroups = []
    found = False
    base_level = None

    for level, haplogroup in haplogroups:
        if haplogroup == haplogroup_name:
            found = True
            base_level = level
            downstream_haplogroups.append(haplogroup)
        elif found and level > base_level:
            downstream_haplogroups.append(haplogroup)
        elif found and level <= base_level:
            break

    return downstream_haplogroups

def extract_ids_for_haplogroups(downstream_haplogroups, input_lines):
    print("正在提取ID...")
    haplogroup_to_ids = {}
    for line in input_lines:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            id, haplogroup = parts
            if haplogroup in downstream_haplogroups:
                if haplogroup not in haplogroup_to_ids:
                    haplogroup_to_ids[haplogroup] = []
                haplogroup_to_ids[haplogroup].append(id)
    return haplogroup_to_ids

def extract_matching_sequences(fasta_path, ids_to_extract):
    print("正在提取符合条件的序列...")
    """从FASTA文件中提取符合条件的序列"""
    sequences = {}
    with open(fasta_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in ids_to_extract:
                sequences[record.id] = record
    return sequences

# ------------------- 新增的序列检查和过滤功能 -------------------
def check_and_filter_sequences(sequences):
    print("正在检查序列质量...")
    """检查序列长度，替换特殊字符，并根据N比例过滤"""
    sequence_lengths = [len(record.seq) for record in sequences.values()]

    # 检查所有序列长度是否一致
    if len(set(sequence_lengths)) > 1:
        print("错误：序列长度不一致，程序终止。")
        exit(1)

    # 检查特殊字符并替换为N
    valid_chars = set("AGCTN-")
    print("正在检查特殊字符并替换为N...")
    for record in sequences.values():
        seq_str = str(record.seq)  # 获取序列的字符串表示
        # 替换非标准碱基字符为N
        new_seq = ''.join([base if base in valid_chars else 'N' for base in seq_str])
        record.seq = Seq(new_seq)  # 使用Bio.Seq模块的Seq类将字符串转换为Seq对象

    # 计算每个序列的N_ratio
    print("正在计算每个序列的N比例...")
    n_ratios = {record.id: record.seq.count('N') / len(record.seq) for record in sequences.values()}

    # 统一设定N比例阈值为0.01
    threshold = 0.01

    # 过滤N比例超过阈值的序列
    filtered_sequences = {id: seq for id, seq in sequences.items() if n_ratios[id] <= threshold}

    # 提示哪些序列被移除
    removed_ids = [id for id in sequences if id not in filtered_sequences]
    if removed_ids:
        print(f"以下序列因N比例过高被移除: {', '.join(removed_ids)}")
    else:
        print("所有序列的N比例都低于阈值，无需移除。")

    return filtered_sequences

def write_fasta_file(output_path, sequences):
    """将提取的序列写入新的FASTA文件"""
    print("正在写入FASTA文件...")
    with open(output_path, 'w') as output_file:
        SeqIO.write(sequences.values(), output_file, "fasta")

# ------------------- 核苷酸差异计算和ρ、σ部分 -------------------
def calculate_nucleotide_difference(seq1, seq2):
    """计算两个序列之间的核苷酸差异"""
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

def calculate_rho_sigma(fasta_file):
    try:
        alignment = AlignIO.read(fasta_file, "fasta")
        if len(alignment) == 0:
            raise ValueError("No sequences found in alignment")

        # 生成序列对
        seq_pairs = [(seq1.seq, seq2.seq) for seq1, seq2 in combinations(alignment, 2)]
        num_comparisons = len(seq_pairs)

        # 如果没有有效的序列对，直接返回None
        if num_comparisons == 0:
            raise ValueError("No valid sequence pairs for comparison")

        # 使用多进程计算序列差异
        with Pool() as pool:
            differences = pool.starmap(calculate_nucleotide_difference, seq_pairs)

        total_differences = sum(differences)
        total_squared_differences = sum(diff ** 2 for diff in differences)

        # 计算rho和sigma
        rho = total_differences / num_comparisons
        sigma_squared = total_squared_differences / (num_comparisons ** 2)
        sigma = math.sqrt(sigma_squared)

        return rho, sigma

    except ValueError as e:
        print(f"Error processing {fasta_file}: {str(e)}")
        return None, None



def remove_bases(seq, positions):
    for pos in sorted(positions, reverse=True):
        seq = seq[:pos-1] + seq[pos:]
    return seq

def write_fasta(fasta_sequences, output_file):
    with open(output_file, 'w') as f:
        for header, seq in fasta_sequences:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(f"{seq[i:i+60]}\n")

def process_fasta(input_file):
    fasta_sequences = list(SeqIO.parse(input_file, 'fasta'))

    output_files = {
        "Complete_sequence.fasta": [],
        "HVS-I_(16051-16400).fasta": [],
        "HVS-II_(68-263).fasta": [],
        "Control_Region.fasta": []
    }
    positions_to_remove = [16519, 16194, 16183, 16182]
    for fasta in fasta_sequences:
        header, sequence = fasta.id, str(fasta.seq)

        seq_1 = remove_bases(sequence, positions_to_remove)
        output_files["Complete_sequence.fasta"].append((header, seq_1))

        seq_2 = seq_1[16051-1:16397]
        output_files["HVS-I_(16051-16400).fasta"].append((header, seq_2))

        seq_3 = sequence[68-1:263]
        output_files["HVS-II_(68-263).fasta"].append((header, seq_3))

        seq_4 = sequence[16024-1:16569] + sequence[0:576]
        output_files["Control_Region.fasta"].append((header, seq_4))

    for output_file, sequences in output_files.items():
        write_fasta(sequences, output_file)

def calculate_for_all_files(output_txt_file, haplogroup_name):
    fasta_files = [
        'Complete_sequence.fasta',
        'HVS-I_(16051-16400).fasta',
        'HVS-II_(68-263).fasta',
        'Control_Region.fasta'
    ]

    results = {}

    for fasta_file in fasta_files:
        print("正在计算，请稍后")
        rho, sigma = calculate_rho_sigma(fasta_file)
        if rho is None or sigma is None:
            print(f"跳过文件 {fasta_file}，因为没有找到有效的序列或比较对")
            continue  # 跳过当前fasta文件

        region_key = fasta_file.split('_')[0].upper()
        results[f'RHO_{region_key}'] = f"{rho:.6f}"
        results[f'RHO_{region_key}_SE'] = f"{sigma:.6f}"

    # 如果没有有效结果，跳过写入文件
    if not results:
        print(f"跳过单倍群 {haplogroup_name}，因为没有计算出有效结果")
        return

    # 使用haplogroup_name作为SAMPLE列
    with open(output_txt_file, 'w') as f:
        f.write("SAMPLE\tRHO_CS\tRHO_CS_SE\tRHO_SYN\tRHO_SYN_SE\t"
                "RHO_HVSI\tRHO_HVSI_SE\tRHO_HVSI_TRANSI\t"
                "RHO_HVSI_TRANSI_SE\tRHO_HVSII\tRHO_HVSII_SE\t"
                "RHO_CR\tRHO_CR_SE\n")
        f.write(f"{haplogroup_name}\t{results.get('RHO_COMPLETE', '0')}\t"
                f"{results.get('RHO_COMPLETE_SE', '0')}\t0\t0\t"
                f"{results.get('RHO_HVS-I', '0')}\t"
                f"{results.get('RHO_HVS-I_SE', '0')}\t0\t0\t"
                f"{results.get('RHO_HVS-II', '0')}\t"
                f"{results.get('RHO_HVS-II_SE', '0')}\t"
                f"{results.get('RHO_CONTROL', '0')}\t"
                f"{results.get('RHO_CONTROL_SE', '0')}\n")

    for fasta_file in fasta_files:
        Path(fasta_file).unlink()


# ------------------- 共祖时间计算部分 -------------------
def calculate_ages(df, rho_col, se_col, multiplier):
    df[f'{rho_col}_AGE'] = df[rho_col] * multiplier
    df[f'{rho_col}_AGE_95LB'] = (df[rho_col] - 1.96 * df[se_col]) * multiplier
    df[f'{rho_col}_AGE_95HB'] = (df[rho_col] + 1.96 * df[se_col]) * multiplier

def process_rho_results(output_txt_file, final_output_file, haplogroup_name, sequence_count):
    df = pd.read_csv(output_txt_file, sep='\t')

    required_columns = [
        'Haplogroup', 'RHO_CS', 'RHO_CS_SE',
        'RHO_SYN', 'RHO_SYN_SE',
        'RHO_HVSI', 'RHO_HVSI_SE',
        'RHO_HVSI_TRANSI', 'RHO_HVSI_TRANSI_SE',
        'RHO_HVSII', 'RHO_HVSII_SE',
        'RHO_CR', 'RHO_CR_SE'
    ]

    for col in required_columns:
        if col not in df.columns:
            df[col] = 0

    # 插入用户输入的单倍群名称
    df['Haplogroup'] = haplogroup_name

    # 插入序列数量
    df['Number'] = sequence_count

    cs_multiplier = 3624
    cs_base = -40.2789
    cs_exp = -0.0263
    df['CS_AGE'] = cs_multiplier * (((np.exp(-(np.exp((df['RHO_CS'] - cs_base) * cs_exp)))) * 0.4794) / 0.4794 * df['RHO_CS'])
    df['CS_AGE_95LB'] = cs_multiplier * (((np.exp(-(np.exp(((df['RHO_CS'] - 1.96 * df['RHO_CS_SE']) - cs_base) * cs_exp)))) * 0.4794) / 0.4794 * (df['RHO_CS'] - 1.96 * df['RHO_CS_SE']))
    df['CS_AGE_95HB'] = cs_multiplier * (((np.exp(-(np.exp(((df['RHO_CS'] + 1.96 * df['RHO_CS_SE']) - cs_base) * cs_exp)))) * 0.4794) / 0.4794 * (df['RHO_CS'] + 1.96 * df['RHO_CS_SE']))

    multipliers = {
        'RHO_SYN': 7872,
        'RHO_HVSI': 16677,
        'RHO_HVSI_TRANSI': 18845,
        'RHO_HVSII': 22388,
        'RHO_CR': 9058
    }

    for rho_col, multiplier in multipliers.items():
        calculate_ages(df, rho_col, f'{rho_col}_SE', multiplier)

    result_columns = [
        'Haplogroup', 'Number', 'CS_AGE', 'CS_AGE_95LB', 'CS_AGE_95HB',
        'RHO_SYN_AGE', 'RHO_SYN_AGE_95LB', 'RHO_SYN_AGE_95HB',
        'RHO_HVSI_AGE', 'RHO_HVSI_AGE_95LB', 'RHO_HVSI_AGE_95HB',
        'RHO_HVSI_TRANSI_AGE', 'RHO_HVSI_TRANSI_AGE_95LB', 'RHO_HVSI_TRANSI_AGE_95HB',
        'RHO_HVSII_AGE', 'RHO_HVSII_AGE_95LB', 'RHO_HVSII_AGE_95HB',
        'RHO_CR_AGE', 'RHO_CR_AGE_95LB', 'RHO_CR_AGE_95HB'
    ]
    result = df[result_columns]
    result.to_csv(final_output_file, sep='\t', index=False, encoding='UTF-8', mode='a', header=False)
    print(result)

# ------------------- 主程序 -------------------
def main():
    # 读取母系单倍群文件
    with open(HAPLOGROUP_FILE_PATH, 'r', encoding='utf-8') as file:
        refined_lines = file.readlines()

    haplogroups = parse_haplogroup_file(refined_lines)

    # 从CSV文件读取单倍群名称
    haplogroup_df = pd.read_csv('C:/Users/victo/Desktop/列表.csv')
    haplogroup_list = haplogroup_df['haplogroup'].tolist()  # 假设CSV文件有一列名为 'haplogroup'

    for haplogroup_name in haplogroup_list:
        print(f"处理单倍群：{haplogroup_name}")

        # 查找该单倍群及其所有下游单倍群
        downstream_haplogroups = find_downstream_haplogroups(haplogroup_name, haplogroups)

        # 读取输入的txt文件
        with open(INPUT_TXT_PATH, 'r', encoding='utf-8') as input_file:
            input_lines = input_file.readlines()

        # 提取下游单倍群对应的ID
        haplogroup_to_ids = extract_ids_for_haplogroups(downstream_haplogroups, input_lines)

        # 收集所有的ID
        all_ids = []
        for ids in haplogroup_to_ids.values():
            all_ids.extend(ids)

        # 从FASTA文件中提取匹配的序列
        matching_sequences = extract_matching_sequences(FASTA_FILE_PATH, all_ids)

        # 进行序列检查和过滤
        filtered_sequences = check_and_filter_sequences(matching_sequences)

        # 输出新的FASTA文件
        global OUTPUT_FASTA_PATH
        OUTPUT_FASTA_PATH = os.path.join(DESKTOP_PATH, f'{haplogroup_name}_sequences.fasta')
        write_fasta_file(OUTPUT_FASTA_PATH, filtered_sequences)

        print(f"已将符合条件的序列输出到: {OUTPUT_FASTA_PATH}")

        # 处理生成的FASTA文件并计算ρ和σ
        process_fasta(OUTPUT_FASTA_PATH)

        # 调用 calculate_for_all_files 时传递 haplogroup_name 参数
        calculate_for_all_files(OUTPUT_TXT_FILE, haplogroup_name)

        # 获取过滤后FASTA文件中的序列数量
        sequence_count = len(filtered_sequences)

        # 处理 output_results.txt 文件并计算共祖时间，添加序列数量和单倍群名称
        process_rho_results(OUTPUT_TXT_FILE, FINAL_OUTPUT_FILE, haplogroup_name, sequence_count)

        # 追加模式将 OUTPUT_TXT_FILE 的内容复制到 ρ结果.txt 文件
        with open(OUTPUT_TXT_FILE, 'r') as source_file:
            with open(RHO_RESULT_FILE, 'a') as target_file:
                shutil.copyfileobj(source_file, target_file)

        print(f"已将 {OUTPUT_TXT_FILE} 的内容追加到 {RHO_RESULT_FILE}")
        # 删除产生的fasta序列
        os.remove(OUTPUT_FASTA_PATH)
        print(f"已删除 {OUTPUT_FASTA_PATH}")

if __name__ == "__main__":
    main()
