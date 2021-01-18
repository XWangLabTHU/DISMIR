# DISMIR
These codes are core codes of DISMIR:

DISMIR_training.py builds and trains a deep learning model to predict the source of individual cell-free DNA reads. The input data files should contain following information:
* the first column : the chromosome where the read is from (e.g. chr1)
* the second column : the position of the read on chromosome (int format)
* the third column of each row in the file should be the sequence of a read with a length of 66bp (A/T/C/G)
* the fourth column of each row in the file should be the methylation state corresponding to the sequence also with a length of 66bp (1:methylated/0:unmethylated) both columns in the form of 'str'
    
DISMIR_value_sample.py predicts the source of each individual read with the trained model from DISMIR_training.py. The input files should have the same format as DISMIR_training.py.

DISMIR_cal_risk.py calculates the risk that the plasma donor is suffering from cancer with the results from DISMIR_value_sample.py

Information about the environment we used:
Python 3.8.3
Tensorflow 2.2.0
Keras 2.3.1


For testing the pipline, we generated several fictitious samples and put them in the directory "train_data" and "test_data". Putting them in the corresponding directory in DISMIR_training.py and DISMIR_value_sample.py, the core codes can be tested.

We also provided scripts to process .BAM files mapped with BS-seeker2. By running the script in "script_for_bam" in order, you can get the .BED file containing switching regions and .read files as the input of the deep learning model. You can process your own data with these scripts and then test the core codes of DISMIR. Here we need samtools to process the .BAM files. Python package pyfaidx and pysam are also required.
