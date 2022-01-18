import numpy as np
import os
import pysam
from pyfaidx import Fasta


bin_length = 500 # bin length of tested regions on the genome
file_dir = '/data/jqli/HCC/DISMIR_gen_data/reads_methy_ratio/' # file dir to store results from training data
bam_dir = '/data/jqli/HCC/DISMIR_gen_data/training_bam/' # file dir of training data bam
hg19_dir = '/data/jqli/hg19/hg19.fa' # the reference genome


def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        return files


if __name__ == '__main__':
    files_0 = file_name(bam_dir)
    genes = Fasta(hg19_dir)
    files = []
    for file in files_0:
        if ('.bai' not in file) & ('.bam' in file):
            files.append(file)

    for file in files:
        bam_temp = pysam.Samfile(bam_dir + file, "rb")
        output = open(file_dir + file + '_' + str(bin_length) + '_m.txt', 'w')
        for i in range(22):
            num = str(i + 1)
            chr_name = 'chr' + num
            length = genes[chr_name][-1].start  # length of current chrom
            region_num = int(np.ceil(length / bin_length))
            for j in range(region_num):
                ratio_H, ratio_L, count = 0, 1, 0
                for read in bam_temp.fetch(chr_name, j * bin_length + 1, (j + 1) * bin_length):
                    unC = read.tags[3][1].count("x")
                    meC = read.tags[3][1].count("X")
                    if (unC + meC) > 3:
                        ratio = meC / (meC + unC)
                        if ratio > ratio_H:
                            ratio_H = ratio
                        if ratio < ratio_L:
                            ratio_L = ratio
                        count = count + 1
                output.write(chr_name + '\t' + str(j * bin_length + 1) + '\t' + str((j + 1) * bin_length) + '\t' + str(
                    ratio_H) + '\t' + str(ratio_L) + '\t' + str(count) + '\n')
        output.close()

