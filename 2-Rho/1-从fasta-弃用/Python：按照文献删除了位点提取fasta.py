from Bio import SeqIO
# 调用函数处理FASTA文件
input_file = r'C:/Users/victo/Desktop/Haplogroup_D.fasta'
def remove_bases(sequence, positions):
    """从序列中移除特定位置的碱基（使用1-based索引）"""
    sequence_list = list(sequence)
    for pos in sorted(positions, reverse=True):
        if 0 < pos <= len(sequence_list):
            sequence_list.pop(pos - 1)  # 转换为0-based索引移除碱基
    return ''.join(sequence_list)

def write_fasta(sequences, output_file):
    """将序列写入到FASTA文件"""
    with open(output_file, 'w') as f:
        for header, sequence in sequences:
            f.write(f">{header}\n")
            f.write(f"{sequence}\n")

def process_fasta(input_file):
    """处理FASTA文件并生成5个新的FASTA文件"""
    # 读取FASTA文件中的所有序列
    fasta_sequences = list(SeqIO.parse(input_file, 'fasta'))

    # 初始化5个输出文件
    output_files = {
        "Complete_sequence.fasta": [],  # 完整序列
        "HVS-I_(16051-16400).fasta": [],  # HVS-I区域
        "HVS-II_(68-263).fasta": [],  # HVS-II区域
        "Control_Region.fasta": [],  # 控制区
        "577_to_16023.fasta": []  # 新的提取区域
    }

    # 需要移除的碱基位置（1-based）
    positions_to_remove = [16519, 16194, 16183, 16182]

    for fasta in fasta_sequences:
        header, sequence = fasta.id, str(fasta.seq)

        # 第一种情况：完整序列移除特定碱基
        seq_1 = remove_bases(sequence, positions_to_remove)
        output_files["Complete_sequence.fasta"].append((header, seq_1))

        # 第二种情况：HVS-I区域（16051-16400）
        seq_2 = seq_1[16051-1:16397]  # 提取HVS-I区间
        output_files["HVS-I_(16051-16400).fasta"].append((header, seq_2))

        # 第三种情况：HVS-II区域（68-263）
        seq_3 = sequence[68-1:263]  # 提取HVS-II区间
        output_files["HVS-II_(68-263).fasta"].append((header, seq_3))

        # 第四种情况：控制区（16024-16569加上1-576）
        seq_4 = sequence[16024-1:16569] + sequence[0:576]  # 提取控制区
        output_files["Control_Region.fasta"].append((header, seq_4))

        # 第五种情况：提取577到16023的序列
        seq_5 = sequence[577-1:16023]  # 提取577-16023区域
        output_files["577_to_16023.fasta"].append((header, seq_5))

    # 写入生成的5个FASTA文件
    for output_file, sequences in output_files.items():
        write_fasta(sequences, output_file)


process_fasta(input_file)
