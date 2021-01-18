import numpy as np
import os
import sys


def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        return files


file_dir = '/data/jqli/HCC/DISMIR_gen_data/reads_methy_ratio/' # file dir of results from training data
root_dir = '/data/jqli/HCC/DISMIR_gen_data/'


if __name__ == '__main__':
    files = file_name(file_dir)
    count = len(open(file_dir + files[0], 'r').readlines())

    normal_P_low = np.ones(count)
    normal_P_num = np.zeros(count)

    tumor_T_low = np.ones(count)
    tumor_T_num = np.zeros(count)

    for file in files:
        if ('CTR' in file) & ('T-' not in file): # keyword for normal plasma samples for training
            data_temp = open(file_dir + file,'r')
            i = 0
            for lines in data_temp:
                line = lines.split()
                if float(line[4])<normal_P_low[i]:
                    normal_P_low[i] = float(line[4])
                normal_P_num[i] = normal_P_num[i]+int(line[5])
                i = i + 1
        if ('HOT' in file) & ('T-' in file): # keyword for cancer tissue samples for training
            data_temp = open(file_dir + file,'r')
            i = 0
            for lines in data_temp:
                line = lines.split()
                if float(line[4])<tumor_T_low[i]:
                    tumor_T_low[i] = float(line[4])
                tumor_T_num[i] = tumor_T_num[i]+int(line[5])
                i = i + 1

    diff = normal_P_low - tumor_T_low
    B_bed = open(root_dir + 'DISMIR.bed','w')
    file = open(file_dir + files[0], 'r')
    i = 0
    for lines in file:
        line = lines.split()
        if (diff[i]>0.3)&(tumor_T_num[i]>25)&(normal_P_num[i]>25):
            B_bed.write(line[0]+'\t'+line[1]+'\t'+line[2]+'\n')
        i = i+1
