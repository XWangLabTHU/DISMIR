# Note that this script is used to process bam files mapped with BS-Seeker2
# In DISMIR the bam files should be preprocessed, for example, use samtools to fetch reads from selected regions
import pysam
import os

bam_dir = '/data/jqli/HCC/bam_process/bam_files/' # where the bam files of selected regions are stored
read_dir = '/data/jqli/HCC/bam_process/read_files/' # where to store the processed read files


def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        return files


files = file_name(bam_dir)

for file in files:
    if ('.bam' in file) & ('log' not in file):
        input = pysam.AlignmentFile(bam_dir + file, 'rb')
        file = file[0:len(file)-4]
        output = open(read_dir + file+'.reads', 'w')
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
            output.write(name+'\t'+start+'\t'+read+'\t'+methylation+'\n')


