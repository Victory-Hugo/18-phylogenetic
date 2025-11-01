from Bio import AlignIO, SeqIO
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path
import math

# 全局变量，定义输入文件和输出文件路径
INPUT_FASTA_FILE = Path('C:/Users/victo/Desktop/D4b_sequences.fasta')
OUTPUT_TXT_FILE = Path('C:/Users/victo/Desktop/output_results.txt')

# 计算核苷酸差异
def calculate_nucleotide_difference(seq1, seq2):
    """计算两个序列之间的核苷酸差异"""
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

def calculate_rho_sigma(fasta_file):
    """计算群体中个体间的平均核苷酸差异 ρ 和 σ"""
    alignment = AlignIO.read(fasta_file, "fasta")

    seq_pairs = [(seq1.seq, seq2.seq) for seq1, seq2 in combinations(alignment, 2)]
    num_comparisons = len(seq_pairs)

    # 并行处理计算核苷酸差异
    with Pool() as pool:
        differences = pool.starmap(calculate_nucleotide_difference, seq_pairs)

    total_differences = sum(differences)
    total_squared_differences = sum(diff ** 2 for diff in differences)

    rho = total_differences / num_comparisons
    sigma_squared = total_squared_differences / (num_comparisons ** 2)
    sigma = math.sqrt(sigma_squared)

    return rho, sigma

def remove_bases(seq, positions):
    """从后往前删除指定位置的碱基"""
    for pos in sorted(positions, reverse=True):
        seq = seq[:pos-1] + seq[pos:]
    return seq

def write_fasta(fasta_sequences, output_file):
    """将处理后的序列写入FASTA文件"""
    with open(output_file, 'w') as f:
        for header, seq in fasta_sequences:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(f"{seq[i:i+60]}\n")

def process_fasta(input_file):
    """处理FASTA文件并生成4个新的FASTA文件"""
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

        # 第一种情况
        seq_1 = remove_bases(sequence, positions_to_remove)
        output_files["Complete_sequence.fasta"].append((header, seq_1))

        # 第二种情况
        seq_2 = seq_1[16051-1:16397]
        output_files["HVS-I_(16051-16400).fasta"].append((header, seq_2))

        # 第三种情况
        seq_3 = sequence[68-1:263]
        output_files["HVS-II_(68-263).fasta"].append((header, seq_3))

        # 第四种情况
        seq_4 = sequence[16024-1:16569] + sequence[0:576]
        output_files["Control_Region.fasta"].append((header, seq_4))

    # 写入生成的4个FASTA文件
    for output_file, sequences in output_files.items():
        write_fasta(sequences, output_file)

def calculate_for_all_files(output_txt_file):
    """对4个新生成的FASTA文件分别进行ρ和σ的计算，并将结果写入txt文件"""
    fasta_files = [
        'Complete_sequence.fasta',
        'HVS-I_(16051-16400).fasta',
        'HVS-II_(68-263).fasta',
        'Control_Region.fasta'
    ]

    results = {}

    for fasta_file in fasta_files:
        rho, sigma = calculate_rho_sigma(fasta_file)
        region_key = fasta_file.split('_')[0].upper()
        results[f'RHO_{region_key}'] = f"{rho:.6f}"
        results[f'RHO_{region_key}_SE'] = f"{sigma:.6f}"

    # 写入结果到txt文件
    with open(output_txt_file, 'w') as f:
        f.write("SAMPLE\tRHO_CS\tRHO_CS_SE\tRHO_SYN\tRHO_SYN_SE\t"
                "RHO_HVSI\tRHO_HVSI_SE\tRHO_HVSI_TRANSI\t"
                "RHO_HVSI_TRANSI_SE\tRHO_HVSII\tRHO_HVSII_SE\t"
                "RHO_CR\tRHO_CR_SE\n")
        f.write(f"1\t{results.get('RHO_COMPLETE', '0')}\t"
                f"{results.get('RHO_COMPLETE_SE', '0')}\t0\t0\t"
                f"{results.get('RHO_HVS-I', '0')}\t"
                f"{results.get('RHO_HVS-I_SE', '0')}\t0\t0\t"
                f"{results.get('RHO_HVS-II', '0')}\t"
                f"{results.get('RHO_HVS-II_SE', '0')}\t"
                f"{results.get('RHO_CONTROL', '0')}\t"
                f"{results.get('RHO_CONTROL_SE', '0')}\n")

    # 删除临时产生的FASTA文件
    for fasta_file in fasta_files:
        Path(fasta_file).unlink()
        print(f"已经删除了临时文件： {fasta_file}，如需保留删掉相应代码")

if __name__ == '__main__':
    # 处理输入FASTA文件并生成新的FASTA文件
    process_fasta(INPUT_FASTA_FILE)

    # 对新生成的FASTA文件计算 ρ 和 σ 并写入结果文件
    calculate_for_all_files(OUTPUT_TXT_FILE)
