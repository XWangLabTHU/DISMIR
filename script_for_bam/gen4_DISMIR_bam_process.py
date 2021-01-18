import pysam
import numpy as np
from collections import Counter
import random

import os
import os.path
import re
import sys
import codecs
import collections


processed_bam_dir = '/data/jqli/HCC/DISMIR_gen_data/processed_bam/' # dir of the processed bam file
reads_dir = '/data/jqli/HCC/DISMIR_gen_data/reads_files/' # dir to store input files of DISMIR deep learning model


def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        return files


def reverse_seq(read, methy):
    read = list(read)
    methy = list(methy)
    methy = methy[1:len(methy)]
    methy.append('.')
    read = "".join(read)
    methy = "".join(methy)
    return read, methy


if __name__ == '__main__':
    files = file_name(processed_bam_dir)
    for file in files:
        if ('.bam' in file) & ('log' not in file):
            input = pysam.AlignmentFile(processed_bam_dir + file, 'rb')
            file = file[0:len(file)-4]
            output = open(reads_dir + file+'.reads', 'w')
            for line in input:
                sig = line.cigarstring
                if 'I' in sig:
                    continue
                name = line.reference_name
                start = str(line.reference_start)
                end = str(line.reference_end)
                read = line.get_tag('XG')
                if 'N' in read:
                    continue
                methy = line.get_tag('XM')
                length = len(read)
                read = read[3:length - 3]
                methylation = ''
                for i in range(len(methy)):
                    if methy[i] == 'X':
                        if i < (len(methy)-1):
                            if (read[i] == 'C') & (read[i+1] == 'G'): # double check
                                methylation = methylation + '1'
                            else:
                                methylation = methylation + '0'
                        elif i == (len(methy)-1):
                            if read[i] == 'C':
                                methylation = methylation + '1'
                            else:
                                methylation = methylation + '0'
                        else:
                            methylation = methylation + '0'
                    else:
                        methylation = methylation+'0'
                if len(methylation) < 71:
                    continue
                else:
                    output.write(name+'\t'+start+'\t'+read[5:71]+'\t'+methylation[5:71]+'\n')

