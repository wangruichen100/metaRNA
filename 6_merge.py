#python 8_merge.py -df gaq_metavirus.xlsx -df_tax nr.xlsx -i virus_r2.fas  -n name -o result.xlsx

import pandas as pd
from Bio import SeqIO
import os
import argparse

# 定义命令行参数解析
parser = argparse.ArgumentParser(description="从不同样本文件中合并并计算RPKM和RPM")

# 定义必要的命令行参数
parser.add_argument('-df', '--datafile', help='包含样本信息的Excel文件路径', required=True)
parser.add_argument('-df_tax', '--taxfile', help='包含注释信息Excel文件路径', required=True)
parser.add_argument('-i', '--input_fasta', help='总病毒序列fas文件', required=True)
parser.add_argument('-n', '--name', help='样本名称列名', required=True)
parser.add_argument('-o', '--output', help='输出文件路径', required=True)

# 解析参数
args = parser.parse_args()

# 读取命令行参数
info_file = args.datafile    # 样本信息Excel文件路径
tax_file = args.taxfile  # 注释信息Excel文件路径
fas_file = args.input_fasta #总病毒序列fas文件
sample_name_column = args.name  # 样本名称的列名
output_file = args.output    # 输出的Excel文件路径

# 读取 sample_id 列表
sample_id = pd.read_excel(info_file)[sample_name_column]

# 初始化结果 DataFrame
result = pd.DataFrame()

df_tax_info = pd.read_excel(tax_file)

# 解析 fasta 文件并建立序列字典
seq_dict = {record.name: record.seq for record in SeqIO.parse(fas_file, "fasta")}
    
# 提取对应的序列，并将其添加到 df_tax_info
df_seq = [seq_dict.get(name, None) for name in df_tax_info["Name"]]
df_tax_info["Seq"] = df_seq

for i in sample_id:
    # 定义文件路径
    df_nr_v_idxstats_path = os.path.join(f"./{i}/bwa/idxstats.txt")
    read_num_path = os.path.join(f"./{i}/1_fastp/{i}_read_num.txt")

    if not os.path.exists(df_nr_v_idxstats_path):
        print(f"跳过...",i)
        continue  # 跳过当前循环
    elif os.path.getsize(df_nr_v_idxstats_path) == 0:
        print(f"跳过...",i)
        continue  # 跳过当前循环

    # 读取 idxstats 文件
    df_nr_v_idxstats = pd.read_table(df_nr_v_idxstats_path, names=["Name", "read_len", "read_num", "non"])

    if not os.path.exists(read_num_path):
        print(f"跳过...",i)
        continue  # 跳过当前循环
    elif os.path.getsize(read_num_path) == 0:
        print(f"跳过...",i)
        continue  # 跳过当前循环


    # 合并 nr_v_info 和 idxstats 数据
    df = df_tax_info.merge(df_nr_v_idxstats, on='Name')
    df["sample"] = [i] * len(df)
    
    # 读取去除rRNA后的总reads数
    read_all_rmrrna = pd.read_csv(read_num_path, sep='\s+')
    read_all_rmrrna['num_seqs'] = read_all_rmrrna['num_seqs'].str.replace(',', '').astype(float)
    df["rmrrna_read"] = [read_all_rmrrna["num_seqs"].sum()] * len(df)

    # 将结果合并到总结果中
    result = pd.concat([result, df])

# 替换所有 0 为 1，防止分母为 0
result = result.replace(0, 1)

# 获取列用于计算 RPKM 和 RPM
contigs_len = result["read_len"].astype(float)
contigs_mapping_reads = result["read_num"].astype(float)
total_mapping_reads = result["rmrrna_read"].astype(float)

# 初始化 RPKM 和 RPM 列表
rpkm_all = []
rpm_all = []

# 计算 RPKM 和 RPM
for x, y, z in zip(contigs_len, contigs_mapping_reads, total_mapping_reads):
    rpkm = (y * 10**9) / (z * x)
    rpm = (y * 10**6) / z
    rpkm_all.append(rpkm)
    rpm_all.append(rpm)

# 将计算结果添加到 result DataFrame
result["rpkm"] = rpkm_all
result["rpm"] = rpm_all

# 保存结果到 Excel 文件
result.to_excel(output_file, index=False)

print(f"结果已保存到 '{output_file}'")
