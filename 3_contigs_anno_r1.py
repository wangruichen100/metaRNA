#python 3_contigs_anno_r1.py -i all_contigs.fas -o nr_v

import os
import re
import argparse
import datetime

# 定义 contigs2anno 函数
def contigs2anno(contigs_path, out_path, db_path, mdb_path):

    # 运行 DIAMOND 和 MEGAN 的系统命令
    print("DIAMOND START: ", datetime.datetime.now())
    os.system(f"diamond blastx --db {db_path} --query {contigs_path} --out {out_path} --id 40 --query-cover 40 --evalue 1e-5 -p 36 --tmpdir /dev/shm --more-sensitive -f 100")
    os.system(f"daa-meganizer -i {out_path}.daa -mdb {mdb_path} --longReads")
    os.system(f"daa2info -i {out_path}.daa -o {out_path}_info.txt -r2c Taxonomy --paths true --list true -mro true")
    
    print("DIAMOND END: ", datetime.datetime.now())

# 主程序：解析命令行参数并处理样本
def main():
    # 定义命令行参数解析
    parser = argparse.ArgumentParser(description='contigs注释')

    # 添加参数
    parser.add_argument('-i', '--contigs_path', type=str, required=True, help='Path to the contigs file')
    parser.add_argument('-o', '--out_path', type=str, required=True, help='Path for the output files')
    
    # 添加具有默认值的参数
    parser.add_argument('--db_path', type=str, default="/mnt/e/ngs/dataset/nr_data/nr_v.dmnd", help='Path to the DIAMOND database (default: /mnt/e/ngs/dataset/nr_data/nr_v.dmnd)')
    parser.add_argument('--mdb_path', type=str, default="/mnt/e/ngs/dataset/megan/megan-map-Feb2022.db", help='Path to the MEGAN mapping database (default: /mnt/e/ngs/dataset/megan/megan-map-Feb2022.db)')

    # 解析参数
    args = parser.parse_args()

    # 调用 contigs2anno 函数
    contigs2anno(args.contigs_path, args.out_path, args.db_path, args.mdb_path)

if __name__ == "__main__":
    main()
