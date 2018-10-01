import sys
import numpy as np

def main():
    args = sys.argv
    # args[1]:threshold number word [2]:number of topic
    #     [3]:cancer_type
    K = int(args[2])
    if(K <= 9):
        topic = '0' + args[2]
    else:
        topic = args[2]

    input_file = 'result/data3_o' + args[1] + '_' + args[3] \
            + '/result_k' + topic + '.txt'
    n4_data = load_data(input_file, K)

    n1_data = to_n1(n4_data, K)
    output_file = 'result/data3_o' + args[1] + '_' + args[3] \
            + '/figure/' + args[2] + '/k' + topic + '.txt'
    write_data(output_file, n1_data, K)

def load_data(input_file, K):
    result = open(input_file,'r')
    lines = result.readlines()
    count = -2
    for line in lines:
        if(count ==  0):
            probability = line.split(' ')
            data = np.zeros([K,len(probability) - 1])
            for i in range(len(probability) - 1):
                data[count,i] = float(probability[i])
        elif((count > 0) and (count < K)):
            probability = line.split(' ')
            for i in range(len(probability) - 1):
                data[count,i] = float(probability[i])
        count += 1
    result.close()
    return data

def to_n1(n4_data, K):
    data = np.zeros([K,96])
    for k in range(K):
        for i in range(96):
            data[k,i] += n4_data[k,i]
    for k in range(K):
        sum = 0
        for i in range(96):
            sum += data[k,i]
        for i in range(96):
            data[k,i] /= sum
    return data

def write_data(output_file, n1_data, K):
    output = open(output_file, 'w')
    output.write('0\n')
    output.write('0\n')
    for k in range(K):
        for i in range(96):
            if(i != 95):
                string = str(n1_data[k,i]) + ' '
            else:
                string = str(n1_data[k,i]) + '\n'
            output.write(string)
    output.close()

if __name__ == '__main__':
    main()
