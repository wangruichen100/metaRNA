#python 4_ex_virus_r1.py -df nr_v.txt -oe nr_v.xlsx -i all_contigs.fas -o virus_r1.fas

import pandas as pd
import datetime
import argparse
import re
import os

def blast_reshape(info_txt, out_excel):
    df = open(info_txt)
    data = []
    for i in df.readlines():
        if not re.match(r"^#", i):
            tax_key = ["Name", "D", "K", "P", "C", "O", "F", "G", "S"]
            tax_info = ["Name", "unknown","unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown"]
            i_ = re.split("[;\t]", i.replace("\n", ""))
            # 格式化一下，删除开头的空格
            for num, j in enumerate(i_):
                i_[num] = re.sub(r"^\s", "", j)
            # 匹配填充
            for num, taxonomy in enumerate(tax_key):
                re_seq = re.compile("^\[%s\].*" % (taxonomy))
                for tax in i_:
                    if re_seq.match(tax):
                        tax_info[num] = tax.replace("[%s] " % (taxonomy), "")
            tax_info[0] = i_[0]
            data.append(tax_info)
            # print(tax_info)
    df_tax = pd.DataFrame(data, columns=["Name", "D", "K", "P", "C", "O", "F", "G", "S"], index=None)
    df_tax.to_excel(out_excel, index=False)
    return df_tax["Name"].to_list()
    
def extract_and_write_sequences(input_fasta, output_fasta, sequence_ids):
    sequences = set(sequence_ids)  # 将序列ID转换为set加快查找速度
    total_sequences = 0

    # 计算输入文件中的序列数量，以便显示进度
    with open(input_fasta) as f:
        for line in f:
            if line.startswith('>'):
                total_sequences += 1

    processed_sequences = 0
    progress_step = max(1, total_sequences // 100)  # 每1%显示一次进度

    with open(input_fasta) as f:
        with open(output_fasta, "w") as output_file:
            write_sequence = False  # 控制是否写入序列
            for line in f:
                if line.startswith(">"):
                    # 提取登录号
                    match = re.match(r'^>(\S+)', line)
                    if match:
                        login_id = match.group(1)
                        # 检查是否在目标序列ID中
                        if login_id in sequences:
                            output_file.write(line)  # 写入ID行
                            write_sequence = True    # 下一行是序列，允许写入
                        else:
                            write_sequence = False
                    # 更新进度
                    processed_sequences += 1
                    if processed_sequences % progress_step == 0:
                        progress = (processed_sequences / total_sequences) * 100
                        print(f"Progress: {progress:.2f}% ({processed_sequences}/{total_sequences} sequences processed)")
                elif write_sequence:
                    output_file.write(line)  # 写入序列行

    print("Sequence extraction completed.")



# 定义命令行参数解析
parser = argparse.ArgumentParser(description="从注释txt提取序列")

# 定义必要的命令行参数
parser.add_argument('-df', '--datafile', help='包含注释信息的txt', required=True)
parser.add_argument('-oe', '--out_excel', help='txt读取为excel', required=True)
parser.add_argument('-i', '--input_fasta', help='输入fasta', required=True)
parser.add_argument('-o', '--output_fasta', help='输出fasta', required=True)

# 解析参数
args = parser.parse_args()

# 读取命令行参数
info_txt = args.datafile
out_excel = args.out_excel
input_fasta = args.input_fasta
output_fasta = args.output_fasta

seq_name = blast_reshape(info_txt, out_excel)
extract_and_write_sequences(input_fasta, output_fasta, seq_name)










