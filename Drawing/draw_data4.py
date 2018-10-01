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
    #[3]:threshold number word [4]:cancer_type
    K = int(args[2])
    topic_list = range(K)
    mutation_list = range(6)

    dictionary_file = 'data/dictionary/data' + args[1] + '_o' + args[3] \
            + '_' + args[4] + '.txt'
    temp_dictionary = load_dictionary(dictionary_file)
    V = len(temp_dictionary)
    dictionary = load_dictionary('data/dictionary.txt')

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
    
    p_ex_proj = np.zeros([K, 1536])
    for v in range(V):
        i = dictionary.index(temp_dictionary[v])
        for k in range(K):
            p_ex_proj[k, i] = p_ex[k, v]

    p_ex_copy = np.zeros([K, 1536])
    pre_labels = [0 for i in range(1536)]
    for i in range(1536):
        fifth = i % 4
        fourth = (i // 4) % 4
        third = ((i // 4) // 4) % 6
        second = (((i // 4) // 4) // 6) % 4
        first = (((i // 4) // 4) // 6) // 4
        j = third * 256 + second * 64 + fifth * 16 + first * 4 + fourth
        pre_labels[j] = make_vocab(first, second, third, fourth, fifth)
        for k in range(K):
            p_ex_copy[k, j] = p_ex_proj[k, i]
    prob = p_ex_copy.copy()

    substitute_labels = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

    for i in range(K):
        if(i in topic_list):
            pre_height = prob[i]
            for j in range(6):
                if(j in mutation_list):
                    fig = plt.figure(figsize=(15,5))
                    left = range(1, 512, 2)
                    height = [0 for l in range(256)]
                    labels = [0 for l in range(256)]
                    l = 0
                    for k in range(j*256, (j+1)*256):
                        height[l] = pre_height[k]
                        labels[l] = pre_labels[k]
                        l += 1
                    title = 'Predicted Signature ' + str(i+1) + ' ' \
                            + substitute_labels[j] + ' in ' + args[4]
                    ax1 = fig.add_subplot(1,1,1)
                    if(j == 0):
                        co = 'r'
                    elif(j == 1):
                        co = 'g'
                    elif(j == 2):
                        co = 'b'
                    elif(j == 3):
                        co = 'c'
                    elif(j == 4):
                        co = 'm'
                    elif(j == 5):
                        co = 'y'
                    ax1.bar(left, height, width = 1.8, \
                            color = co, align='center')
                    max_height = max(height)
                    height_lim = math.ceil(max_height*100)/100
                    ax1.set_ylim(0,height_lim)
                    ax1.set_xlim(0,514)
                    ax1.set_xticks(left)
                    ax1.set_xlabel('mutation x', fontsize=9)
                    ax1.set_xticklabels(labels)
                    for tick in ax1.get_xticklabels():
                        tick.set_rotation(90)
                    ax1.tick_params(labelsize = 3)
                    ax1.set_ylabel('p (mutation = x)', fontsize=9)
                    ax1.set_title(title, fontsize=9)
                    fig.tight_layout()
                    name = 'result/data'+ args[1] + '_o' + args[3] + '_' \
                            + args[4] + '/figure/' + args[2] \
                            + '/detail_context_' + str(i+1) + '_' \
                            + substitute_labels[j][0] + 'to' \
                            + substitute_labels[j][2] + '.png'
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

def make_vocab(first, second, third, fourth, fifth):
    vocab = base(first)
    vocab += base(second)
    if(third == 0):
        vocab += '(C>A)'
    elif(third == 1):
        vocab += '(C>G)'
    elif(third == 2):
        vocab += '(C>T)'
    elif(third == 3):
        vocab += '(T>A)'
    elif(third == 4):
        vocab += '(T>C)'
    elif(third == 5):
        vocab += '(T>G)'
    vocab += base(fourth)
    vocab += base(fifth)
    return vocab

def base(num):
    if(num == 0):
        return('A')
    elif(num == 1):
        return('C')
    elif(num == 2):
        return('G')
    elif(num == 3):
        return('T')
    else:
        return('?')

if __name__ == '__main__':
    main()
