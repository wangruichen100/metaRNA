#python 7_quantify.py -df info.xlsx -n name -1 r1 -2 r2 -i all_virus.fas

from Bio import SeqIO
import pandas as pd
import datetime
import argparse
import re
import os

def extract_and_write_sequences(input_fasta, output_fasta, sample):
    name_contains = sample + "@"
    with open(input_fasta) as f:
        with open(output_fasta, "w") as output_file:
            write_sequence = False  # 控制是否写入序列
            for line in f:
                if line.startswith(">"):
                    # 提取登录号
                    match = re.match(r'^>(\S+)', line)
                    if match:
                        login_id = match.group(1)
                        # 检查登录号中是否包含传入的子字符串
                        if name_contains in login_id:
                            output_file.write(line)  # 写入ID行
                            write_sequence = True    # 下一行是序列，允许写入
                        else:
                            write_sequence = False
                elif write_sequence:
                    output_file.write(line)  # 写入序列行

def quantify(sample, sample1, sample2, all_virus):
    print("{sample} START: ".format(sample = sample), datetime.datetime.now())

    os.system("seqkit stat ./{sample}/1_fastp/*clean.fastq.gz > ./{sample}/1_fastp/{sample}_read_num.txt".format(sample= sample))

    # bwa
    print("bwa START: ", datetime.datetime.now())
    os.system("mkdir ./{sample}/bwa".format(sample= sample))

    extract_and_write_sequences(all_virus,
                                                  "./{sample}/bwa/{sample}_virus.fas".format(sample = sample),
                                                   sample)

    os.system("bwa index ./{sample}/bwa/{sample}_virus.fas".format(sample = sample))
    os.system("bwa mem -t 36 ./{sample}/bwa/{sample}_virus.fas ./{sample}/1_fastp/{sample}_R1_clean.fastq.gz ./{sample}/1_fastp/{sample}_R2_clean.fastq.gz|\
               samtools view -bS -@ 36|\
               samtools sort -@ 36 -o ./{sample}/bwa/bwa_sort.bam".format(sample = sample, sample1 = sample1,sample2 = sample2))
    os.system("samtools index -@ 36 ./{sample}/bwa/bwa_sort.bam".format(sample = sample))
    os.system("samtools flagstat -@ 36 ./{sample}/bwa/bwa_sort.bam > ./{sample}/bwa/flagstat.txt".format(sample = sample))
    os.system("samtools depth ./{sample}/bwa/bwa_sort.bam > ./{sample}/bwa/depth.txt".format(sample = sample))
    os.system("samtools idxstats ./{sample}/bwa/bwa_sort.bam > ./{sample}/bwa/idxstats.txt".format(sample = sample))

    print("{sample} END: ".format(sample = sample), datetime.datetime.now())

parser = argparse.ArgumentParser(description='NGS_VIRUS')
parser.add_argument('-df', help='输入文件路径')
parser.add_argument('-n','--name', help='sample名称')
parser.add_argument('-1','--r1', help='双端测序数据R1路径')
parser.add_argument('-2','--r2', help='双端测序数据R2路径')
parser.add_argument('-i','--input_fas', help='all_virus.fas')

args = parser.parse_args()

info_path = args.df
sample_name_column_name = args.name
sample_R1_column_name = args.r1
sample_R2column_name = args.r2
all_virus = args.input_fas

data_info = pd.read_excel(info_path)

sample_name = data_info[sample_name_column_name].to_list()
sample_R1_list = data_info[sample_R1_column_name].to_list()
sample_R2list = data_info[sample_R2column_name].to_list()

for x,y,z in zip(sample_name,sample_R1_list,sample_R2list):
    quantify(x,y,z,all_virus)
