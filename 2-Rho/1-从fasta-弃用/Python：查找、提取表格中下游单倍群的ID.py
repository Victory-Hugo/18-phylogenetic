import os
from Bio import SeqIO  # 引入biopython的SeqIO模块

# 定义桌面路径
DESKTOP_PATH = os.path.join(os.path.expanduser("~"), 'Desktop')
HAPLOGROUP_FILE_PATH = os.path.join('F:/OneDrive/文档（科研）/脚本/我的科研脚本/Python/母系专用/线粒体单倍群phylotree(version17)2024年8月19日.txt')
INPUT_TXT_PATH = r'C:/Users/victo/Desktop/新建 Text Document.txt'
FASTA_FILE_PATH = r'C:/Users/victo/Desktop/Haplogroup_D.fasta'  # fasta文件路径

def parse_haplogroup_file(lines):
    haplogroups = []
    for line in lines:
        stripped_line = line.strip()
        if stripped_line:
            level = line.count('\t')
            haplogroups.append((level, stripped_line))
    return haplogroups

def find_downstream_haplogroups(haplogroup_name, haplogroups):
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
    """从FASTA文件中提取符合条件的序列"""
    sequences = {}
    with open(fasta_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in ids_to_extract:
                sequences[record.id] = record
    return sequences

def write_fasta_file(output_path, sequences):
    """将提取的序列写入新的FASTA文件"""
    with open(output_path, 'w') as output_file:
        SeqIO.write(sequences.values(), output_file, "fasta")

def main():
    # 读取母系单倍群文件
    with open(HAPLOGROUP_FILE_PATH, 'r', encoding='utf-8') as file:
        refined_lines = file.readlines()

    haplogroups = parse_haplogroup_file(refined_lines)

    # 交互式要求用户输入单倍群名称
    haplogroup_name = input("请输入单倍群名称：")

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

    # 输出新的FASTA文件
    output_fasta_path = os.path.join(DESKTOP_PATH, f'{haplogroup_name}_sequences.fasta')
    write_fasta_file(output_fasta_path, matching_sequences)

    print(f"已将符合条件的序列输出到: {output_fasta_path}")

if __name__ == "__main__":
    main()
