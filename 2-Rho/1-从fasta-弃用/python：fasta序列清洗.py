from Bio import SeqIO

# 输入文件路径
input_file = r"C:/Users/victo/Desktop/古代现代对齐.fasta"
output_file = r"C:/Users/victo/Desktop/Wash.fasta"

# 让用户输入N含量的最大比例
max_n_content = float(input("请输入序列中N值含量的最大比例(0~1): 建议输入0.1"))

# 创建一个字典用于存储符合条件的唯一序列ID
unique_records = {}
# 用于存储被删除的序列及原因
removed_records = []

# 读取FASTA文件，使用SeqIO解析
with open(input_file, "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        # 检查序列长度
        if len(record.seq) < 16000:
            removed_records.append((record.id, "序列长度小于 16000"))
            continue
        
        # 计算序列中'N'的比例
        n_count = record.seq.upper().count('N')
        n_content = n_count / len(record.seq)
        
        # 检查N含量是否高于用户输入的最大比例
        if n_content > max_n_content:
            removed_records.append((record.id, f"N 含量 {n_content:.2%} 高于设定的 {max_n_content:.2%}"))
            continue

        # 如果ID不在字典中，添加到字典
        if record.id not in unique_records:
            unique_records[record.id] = record

# 写入新的FASTA文件，仅包含符合条件的唯一序列
with open(output_file, "w") as output_fasta:
    SeqIO.write(unique_records.values(), output_fasta, "fasta")

print(f"去重并过滤后的序列已保存到 {output_file}")

# 打印出被删除的序列的ID和原因
if removed_records:
    print("\n被删除的序列及原因：")
    for record_id, reason in removed_records:
        print(f"ID: {record_id}, 删除原因: {reason}")
else:
    print("没有任何序列被删除。")
