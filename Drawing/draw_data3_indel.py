import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math

def main():
    args = sys.argv
    #args[1]:data_type [2]:number_of_topic
    #    [3]:threshold number word [4]:cancer_type
    K = int(args[2])
    
    dictionary_file_indel = 'data/dictionary/data3_o' + args[3] + '_' \
            + args[4] + '_indel.txt'
    temp_dictionary_indel = load_dictionary(dictionary_file_indel)
    dictionary_indel = load_dictionary('data/dictionary_indel.txt')

    Sub = 96
    V = Sub + len(temp_dictionary_indel)

    if(K <= 9):
        topic = '0' + args[2]
    else:
        topic = args[2]

    input_file = 'result/data' + args[1] + '_o' + args[3] \
               + '_' + args[4] + '/result_k' + topic + '.txt'

    data_ex = open(input_file).readlines()
        
    count_ex = 0
    p_ex = np.zeros([K,V])
    for line in data_ex:
        if((count_ex > 1) and (count_ex < K + 2)):
            words = line.split()
            for signature in range(V):
                p_ex[count_ex-2,signature] = float(words[signature])
        count_ex += 1
    p_ex_copy = np.zeros([K,V-Sub])
    for k in range(K):
        for i in range(V-Sub):
            p_ex_copy[k,i] = p_ex[k, Sub+i]
    indel_pre = p_ex_copy.copy()

    print(temp_dictionary_indel)
    print(dictionary_indel)

    indel = np.zeros([K, 20])
    for k in range(K):
        for v in range(V-Sub):
            indel[k,dictionary_indel.index(temp_dictionary_indel[v])] \
                = indel_pre[k,v]

    labels = [0 for i in range(20)]
    labels[0] = 'A[big]T'
    labels[1] = 'C[big]G'
    labels[2] = 'G[big]C'
    labels[3] = 'T[big]A'
    labels[4] = 'A[big]C'
    labels[5] = 'C[big]T'
    labels[6] = 'T[big]G'
    labels[7] = 'G[big]A'
    labels[8] = 'A[big]A'
    labels[9] = 'C[big]C'
    labels[10] = 'A[small]T'
    labels[11] = 'C[small]G'
    labels[12] = 'G[small]C'
    labels[13] = 'T[small]A'
    labels[14] = 'A[small]C'
    labels[15] = 'C[small]T'
    labels[16] = 'T[small]G'
    labels[17] = 'G[small]A'
    labels[18] = 'A[small]A'
    labels[19] = 'C[small]C'
    colorlist = [0 for i in range(20)]
    for i in range(20):
        if(i < 10):
            colorlist[i] = 'r'
        else:
            colorlist[i] = 'b'

    for i in range(K):
        height = indel[i]
        fig = plt.figure(figsize=(4,2))
        left = range(1,21,1)
        title = 'Predicted Signature ' + str(i+1) + ' (indel) in ' + args[4]
        ax1 = fig.add_subplot(1,1,1)
        ax1.bar(left, height, width = 1, color = colorlist, align='center')
        max_height = max(height)
        height_lim = math.ceil(max_height*10)/10
        ax1.set_ylim(0, height_lim)
        ax1.set_xlim(0, 21)
        ax1.set_xticks(left)
        ax1.set_xlabel('indel', fontsize=6)
        ax1.set_xticklabels(labels)
        for tick in ax1.get_xticklabels():
            tick.set_rotation(90)
        ax1.tick_params(labelsize = 4)
        ax1.set_ylabel('p (indel = x)', fontsize=6)
        ax1.set_title(title, fontsize=6)
        fig.tight_layout()
        name = 'result/data'+ args[1] + '_o' + args[3] + '_' + args[4] \
                + '/figure/' + args[2] + '/indel_' + str(i+1) + '.png'
        fig.savefig(name, dpi=200)
        plt.close(1)

def load_dictionary(dictionary_file):
    file = open(dictionary_file, 'r')
    lines = file.readlines()
    dictionary = list()
    for line in lines:
        text = line[:-1]
        dictionary.append(text)
    return dictionary

if __name__ == '__main__':
    main()
