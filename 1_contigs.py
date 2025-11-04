#python 1_contigs.py -df info.xlsx -n name -1 r1 -2 r2

from Bio import SeqIO
import pandas as pd
import datetime
import argparse
import re
import os

def blast_reshape(sample):
    df = open("./{sample}/3_blast/nr_v_info.txt".format(sample = sample))
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
    df_tax.to_excel("./{sample}/3_blast/nr_v_info.xlsx".format(sample = sample), index=False)
    return df_tax["Name"].to_list()
    
def extract_and_write_sequences(input_fasta, output_fasta, sequence_ids):
    sequences = {}
    with open(input_fasta) as f:
        with open(output_fasta, "w") as output_file:
            for record in SeqIO.parse(f, "fasta"):
                # 使用正则表达式从记录的ID中提取登录号
                match = re.match(r'^(\S+)', record.id)
                if match:
                    login_id = match.group(1)
                    if login_id in sequence_ids:
                        sequences[record.id] = record.seq
                        output_file.write(f">{record.id}\n{record.seq}\n")

def raw2quantify(sample, sample1, sample2):
    print("{sample} START: ".format(sample = sample), datetime.datetime.now())
    os.system("mkdir {sample}".format(sample = sample))
    
    # fastp
    print("fastp START: ", datetime.datetime.now())
    os.system(f"mkdir ./{sample}/1_fastp")
    os.system(f"fastp -i {sample1} -o ./{sample}/1_fastp/{sample}_R1_clean.fastq.gz \
    -I {sample2} -O ./{sample}/1_fastp/{sample}_R2_clean.fastq.gz \
    -j ./{sample}/1_fastp/{sample}.json \
    -h ./{sample}/1_fastp/{sample}.html \
    -w 16")
    # megahit
    print("megahit START: ", datetime.datetime.now())
    os.system(f"megahit -1  ./{sample}/1_fastp/{sample}_R1_clean.fastq.gz \
    -2  ./{sample}/1_fastp/{sample}_R2_clean.fastq.gz \
    -o ./{sample}/2_contig --tmp-dir /tmp -t 36")
    # quast

    print("{sample} END: ".format(sample = sample), datetime.datetime.now())

parser = argparse.ArgumentParser(description='NGS_VIRUS')
parser.add_argument('-df', help='输入文件路径')
parser.add_argument('-n','--name', help='sample名称')
parser.add_argument('-1','--r1', help='双端测序数据R1路径')
parser.add_argument('-2','--r2', help='双端测序数据R2路径')

args = parser.parse_args()

info_path = args.df
sample_name_column_name = args.name
sample_R1_column_name = args.r1
sample_R2column_name = args.r2

data_info = pd.read_excel(info_path)

sample_name = data_info[sample_name_column_name].to_list()
sample_R1_list = data_info[sample_R1_column_name].to_list()
sample_R2list = data_info[sample_R2column_name].to_list()

for x,y,z in zip(sample_name,sample_R1_list,sample_R2list):
    raw2quantify(x,y,z)
