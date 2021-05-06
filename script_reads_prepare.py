# This script is used to process read file into certain length, in DISMIR, the input reads are with the length of 66
import os


def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        return files


read_dir = '/data/jqli/HCC/bam_process/read_files/' # where the read files generated from bam files are stored
read_same_length_dir = '/data/jqli/HCC/bam_process/read_same_length/' # where to store read files with the same length

files = file_name(read_dir)  # 载入文件夹内所有文件名称
for file in files:
    input = open(read_dir + file, 'r')
    output = open(read_same_length_dir + file, 'w')
    for item in input:
        item = item.split()
        if len(item[3]) < 71:
            continue
        else:
            output.write(item[0] + ' ' + item[1] + ' ' + item[2][5:71] + ' ' + item[3][5:71] + '\n')
    input.close()
    output.close()