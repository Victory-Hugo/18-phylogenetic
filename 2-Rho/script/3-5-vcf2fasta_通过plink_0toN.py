import os
import sys
import shutil
from pathlib import Path


def replace_zeros_with_n_inplace(file_path):
    """
    将FASTA文件中的'0'替换为'N'，原地修改文件
    源文件备份为.bak文件
    """
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()

        # 用于存储修改后的内容
        modified_lines = []
        has_zeros = False

        for line in lines:
            line = line.rstrip('\n\r')

            if line.startswith(">"):  # 这是标题行，不做修改
                modified_lines.append(line)
            else:
                # 替换序列中的 '0' 为 'N'
                if '0' in line:
                    has_zeros = True
                modified_sequence = line.replace('0', 'N')
                modified_lines.append(modified_sequence)

        if not has_zeros:
            print(f"✓ {file_path}: 无需修改（未发现'0'字符）")
            return False

        # 备份原文件
        backup_path = file_path + '.bak'
        shutil.copy2(file_path, backup_path)

        # 将修改后的内容写入原文件
        with open(file_path, 'w', encoding='utf-8') as output_file:
            for modified_line in modified_lines:
                output_file.write(modified_line + "\n")

        print(f"✓ {file_path}: 已修改，备份为 {backup_path}")
        return True

    except Exception as e:
        print(f"✗ 处理文件失败 {file_path}: {str(e)}", file=sys.stderr)
        return False


def process_fasta_directory(directory_path):
    """
    处理指定目录中的所有FASTA文件
    """
    dir_path = Path(directory_path)

    if not dir_path.exists():
        print(f"✗ 目录不存在: {directory_path}", file=sys.stderr)
        sys.exit(1)

    if not dir_path.is_dir():
        print(f"✗ 指定路径不是目录: {directory_path}", file=sys.stderr)
        sys.exit(1)

    # 查找所有FASTA文件（.fasta, .fa, .fna等）
    fasta_extensions = ['.fasta', '.fa', '.fna', '.ffn', '.frn']
    fasta_files = []

    for ext in fasta_extensions:
        fasta_files.extend(dir_path.glob(f'*{ext}'))

    if not fasta_files:
        print(f"⚠ 在 {directory_path} 中未找到FASTA文件")
        return

    print(f"开始处理 {len(fasta_files)} 个FASTA文件...\n")

    success_count = 0
    for fasta_file in sorted(fasta_files):
        if replace_zeros_with_n_inplace(str(fasta_file)):
            success_count += 1

    print(f"\n处理完成！共修改 {success_count}/{len(fasta_files)} 个文件")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("用法: python3 script_name.py <fasta_文件夹路径>")
        print("例如: python3 script_name.py /path/to/fasta/directory")
        sys.exit(1)

    directory = sys.argv[1]
    process_fasta_directory(directory)
