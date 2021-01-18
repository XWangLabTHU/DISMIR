import sys
import os
import numpy as np

root_dir = '/data/jqli/HCC/DISMIR_gen_data/'
bed_name = root_dir + 'DISMIR.bed'
bam_dir = '/data/jqli/cell-free/whole_genome/bam_sorted/' # dir of the original mapped bam file
dest_dir = '/data/jqli/HCC/DISMIR_gen_data/processed_bam/' # dir to store the processed bam file


def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        return files


files = file_name(bam_dir)
for file in files:
    if ('.bam' in file) & ('.bai' not in file) & ('.log' not in file):
        os.system('samtools view -b  -h ' + bam_dir + file + ' -L "' + bed_name + '" >' + dest_dir + file + '')