import numpy as np
import random
import os
import os.path
import re

from keras.models import Sequential
from keras import layers
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras import optimizers
import tensorflow as tf


def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        return files


# transform sequence into number for storage (A/T/C/G to 0/1/2/3, methylated C to 4)
def lstm_seq(seq, methy):
    i = 0
    lstmseq = np.zeros((len(seq), 66), dtype='int')
    while i < len(seq):
        tmp = seq[i]
        j = 0
        while j < len(tmp):
            if tmp[j] == 'A':
                lstmseq[i, j] = 0
            elif tmp[j] == 'T':
                lstmseq[i, j] = 1
            elif tmp[j] == 'C':
                lstmseq[i, j] = 2
            else:
                lstmseq[i, j] = 3
            if int(methy[i][j]) == 1:
                lstmseq[i, j] = 4
            j = j + 1
        i = i + 1
    return lstmseq


# transform sequence into one-hot code (0/1/2/3 to one-hot) and add methylation state channel
def conv_onehot(seq):
    module = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 1, 0, 1]])
    onehot = np.zeros((len(seq), 66, 5), dtype='int')
    for i in range(len(seq)):
        tmp = seq[i]
        tmp_onehot = np.zeros((66, 5), dtype='int')
        for j in range(len(tmp)):
            if tmp[j] == 0:
                tmp_onehot[j] = module[0]
            elif tmp[j] == 1:
                tmp_onehot[j] = module[1]
            elif tmp[j] == 2:
                tmp_onehot[j] = module[2]
            elif tmp[j] == 3:
                tmp_onehot[j] = module[3]
            else:
                tmp_onehot[j] = module[4]
        onehot[i] = tmp_onehot
    return onehot


# deep learning model
def DISMIR_deep():
    model = Sequential()
    model.add(layers.Convolution1D(input_shape=(66, 5),
                                   nb_filter=100,
                                   filter_length=10,
                                   border_mode="same",
                                   activation="relu",
                                   subsample_length=1))
    model.add(layers.MaxPooling1D(pool_length=2, stride=2))
    model.add(layers.Dropout(0.2))
    model.add(layers.Bidirectional(layers.LSTM(33, return_sequences=True)))
    model.add(layers.Convolution1D(input_shape=(33, 132),
                                   nb_filter=100,
                                   filter_length=3,
                                   border_mode="same",
                                   activation="relu",
                                   subsample_length=1))
    model.add(layers.MaxPooling1D(pool_length=2, stride=2))
    model.add(layers.Dropout(0.2))
    model.add(layers.Flatten())
    model.add(layers.Dense(750, activation='relu', kernel_regularizer=None, bias_regularizer=None))
    model.add(layers.Dropout(0.2))
    model.add(layers.Dense(300, activation='relu', kernel_regularizer=None, bias_regularizer=None))
    model.add(layers.Dense(1, activation='sigmoid', kernel_regularizer=None, bias_regularizer=None))
    return model


train_dir = '/data/jqli/HCC/12_22_test_program/train_dir/' # directory where the model is saved
file_dir = '/data/jqli/HCC/4_1_one_hot_11/person_5_71/' # directory where the sample to test is saved, format same as training data
store_dir = '/data/jqli/HCC/12_22_test_program/store_dir/' # directory to save predicted d-scores of reads


if __name__ == '__main__':
    gpus = tf.config.experimental.list_physical_devices(device_type='GPU')
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
    model = DISMIR_deep()
    model.load_weights(train_dir + 'weight.h5')

    files = file_name(file_dir)
    for file in files:
        input = open(file_dir + file,'r')
        seq = []
        methy = []
        for item in input:
            item = item.split()
            cpg = 0
            for i in range(len(item[2]) - 1):
                if (item[2][i] == 'C') & (item[2][i + 1] == 'G'):
                    cpg = cpg + 1
            if cpg > 2:
                seq.append(item[2])
                methy.append(item[3])
        input.close()
        seq_lstm = lstm_seq(seq, methy)
        seq_3one_hot = conv_onehot(seq_lstm)
        result = model.predict(seq_3one_hot, verbose=0)
        np.savetxt(store_dir + 'result_' + file + '.txt', result)








