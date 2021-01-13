import numpy as np
import os


def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        return files


if __name__ == '__main__':
    store_dir = '/data/jqli/HCC/12_22_test_program/store_dir/' # directory where the predicted d-scores are stored
    output = open('/data/jqli/HCC/12_22_test_program/value_result.txt', 'w') # file to store predicted risk
    files = file_name(store_dir)

    # cal risk by maximize posterior probability
    gaps = np.linspace(0, 1, 1001)
    score = np.vstack((gaps, 1 - gaps))
    score = score.T
    for file in files:
        likelihood_1 = np.loadtxt(store_dir + file)
        likelihood_2 = 1 - likelihood_1
        likelihood = np.vstack((likelihood_1, likelihood_2))
        val = np.log10(np.dot(score, likelihood))
        sum = np.sum(val, axis=1)
        result = gaps[np.argmax(sum)]
        output.write(file + '\t' + str(result) + '\n')
    output.close()

