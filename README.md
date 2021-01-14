# DISMIR
These codes are core codes of DISMIR

DISMIR_training.py builds and trains a deep learning model to predict the source of individual cell-free DNA reads. The input data files should contain following information:
the first column : the chromosome where the read is from (e.g. chr1)
the second column : the position of the read on chromosome (int format)
the third column of each row in the file should be the sequence of a read with a length of 66bp (A/T/C/G)
the fourth column of each row in the file should be the methylation state corresponding to the sequence also with a length of 66bp (1:methylated/0:unmethylated) both columns in the form of 'str'
    
DISMIR_value_sample.py predicts the source of each individual read with the trained model from DISMIR_training.py. The input files should have the same format as DISMIR_training.py.

DISMIR_cal_risk.py calculates the risk that the plasma donor is suffering from cancer with the results from DISMIR_value_sample.py
