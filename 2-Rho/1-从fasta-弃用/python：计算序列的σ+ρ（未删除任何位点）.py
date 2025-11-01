from Bio import AlignIO
from itertools import combinations
from multiprocessing import Pool
import math

def calculate_nucleotide_difference(seq1, seq2):
    """
    计算两个序列之间的核苷酸差异
    """
    differences = sum(1 for a, b in zip(seq1, seq2) if a != b)
    return differences

def calculate_rho_sigma(fasta_file):
    """
    计算群体中个体间的平均核苷酸差异 ρ 和 σ
    """
    # 读取FASTA文件并进行序列比对
    alignment = AlignIO.read(fasta_file, "fasta")
    
    # 初始化变量
    total_differences = 0
    total_squared_differences = 0
    
    # 准备所有配对组合
    seq_pairs = [(seq1.seq, seq2.seq) for seq1, seq2 in combinations(alignment, 2)]
    num_comparisons = len(seq_pairs)
    
    # 使用多进程池进行并行计算
    with Pool() as pool:
        differences = pool.starmap(calculate_nucleotide_difference, seq_pairs)
    
    # 计算ρ和σ
    for difference in differences:
        total_differences += difference  # n_i
        total_squared_differences += difference ** 2  # n_i^2
    
    # 根据文献公式计算ρ和σ²
    rho = total_differences / num_comparisons
    sigma_squared = total_squared_differences / (num_comparisons ** 2)
    sigma = math.sqrt(sigma_squared)
    
    return rho, sigma

if __name__ == '__main__':
    # 指定FASTA文件路径
    fasta_file = r'C:/Users/victo/Desktop/output_2.fasta'

    # 计算平均核苷酸差异 ρ 和 σ
    rho, sigma = calculate_rho_sigma(fasta_file)
    print(f"平均核苷酸差异 ρ: {rho}")
    print(f"σ值: {sigma}")
