#python 2_merge_contigs.py -df info.xlsx -n name -o all_contigs.fas

import os
import argparse
import pandas as pd

# 定义参数解析
parser = argparse.ArgumentParser(description='合并contigs文件')
parser.add_argument('-df', help='输入包含样本信息的Excel文件路径', required=True)
parser.add_argument('-n', '--name', help='样本名称所在的列名', required=True)
parser.add_argument('-o', '--out', help='合并后的输出fas文件名', required=True)

# 解析命令行参数
args = parser.parse_args()

# 参数变量
info_path = args.df
sample_name_column_name = args.name
output_file = args.out

# 读取样本信息的Excel文件
try:
    sample_info = pd.read_excel(info_path)
except FileNotFoundError:
    print(f"Excel文件 '{info_path}' 不存在，请检查文件路径。")
    exit()

# 检查是否提供的列名存在
if sample_name_column_name not in sample_info.columns:
    print(f"Excel文件中未找到列 '{sample_name_column_name}'，请检查列名。")
    exit()

# 提取样本名称列表
sample_ids = sample_info[sample_name_column_name].tolist()

# 打开输出文件准备合并内容
with open(output_file, 'w') as outfile:
    for sample_id in sample_ids:
        print(sample_id)
        # 构造contigs.fas文件的路径
        contigs_file_path = os.path.join(f"./{sample_id}/2_contig/final.contigs.fa")
        
        # 检查文件是否存在，如果不存在则跳过
        if not os.path.exists(contigs_file_path):
            print(f"样本 {sample_id} 的 contigs 文件 '{contigs_file_path}' 不存在，跳过该样本。")
            continue
        
        # 读取该文件的内容并合并
        with open(contigs_file_path, 'r') as infile:
            for line in infile:
                # 如果是序列名称行，前面加上sample@
                if line.startswith('>'):
                    line = f">{sample_id}@{line[1:]}"
                outfile.write(line)

print(f"所有有效的 contigs 文件已合并到 '{output_file}'")



