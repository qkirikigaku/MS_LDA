import sys
import numpy as np


def main():
    args = sys.argv
    # args[1]:threshold number word [2]:number of topic
    #    [3]:cancer_type
    K = int(args[2])
    if(K <= 9):
        topic = '0' + args[2]
    else:
        topic = args[2]

    input_file = 'result/data2_o' + args[1] + '_' + args[3] \
        + '/result_k' + topic + '.txt'
    n2_data = load_data(input_file, K)

    dictionary_file = 'data/dictionary/data2_o' + args[1] + '_' \
        + args[3] + '.txt'
    temp_dictionary = load_dictionary(dictionary_file)
    dictionary = load_dictionary('data/dictionary.txt')

    n1_data = to_n1(n2_data, dictionary, temp_dictionary, K)
    output_file = 'result/data2_o' + args[1] + '_' + args[3] \
        + '/figure/' + args[2] + '/k' + topic + '.txt'
    write_data(output_file, n1_data, K)


def load_data(input_file, K):
    result = open(input_file, 'r')
    lines = result.readlines()
    count = -2
    for line in lines:
        if(count == 0):
            probability = line.split(' ')
            data = np.zeros([K, len(probability) - 1])
            for i in range(len(probability) - 1):
                data[count, i] = float(probability[i])
        elif((count > 0) and (count < K)):
            probability = line.split(' ')
            for i in range(len(probability) - 1):
                data[count, i] = float(probability[i])
        count += 1
    result.close()
    return data


def load_dictionary(dictionary_file):
    dictionary = list()
    result = open(dictionary_file, 'r')
    lines = result.readlines()
    for line in lines:
        text = line[:-1]
        dictionary.append(text)
    return dictionary


def to_n1(n2_data, dictionary, temp_dictionary, K):
    V = len(temp_dictionary)
    index = [0 for i in range(V)]
    for k in range(V):
        i = dictionary.index(temp_dictionary[k])
        forward = (i % 384) // 96
        substitution = ((i % 384) % 96) // 16
        backward = (((i % 384) % 96) % 16) // 4
        index[k] = 24 * forward + 4 * substitution + backward
    data = np.zeros([K, 96])
    for k in range(K):
        for i in range(V):
            data[k, index[i]] += n2_data[k, i]
    return data


def write_data(output_file, n1_data, K):
    output = open(output_file, 'w')
    output.write('0\n')
    output.write('0\n')
    for k in range(K):
        for i in range(96):
            if(i != 95):
                string = str(n1_data[k, i]) + ' '
            else:
                string = str(n1_data[k, i]) + '\n'
            output.write(string)
    output.close()


if __name__ == '__main__':
    main()
