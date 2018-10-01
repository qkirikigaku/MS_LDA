import re
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import math


def main():
    K = 30
    co = open('data/signature_probability.txt')

    data_co = co.readlines()
    count_co = 0
    p_co = np.zeros([30, 96])
    for line in data_co:
        if(count_co != 0):
            words = line.split()
            for signature in range(33):
                if(signature >= 3):
                    p_co[signature-3, count_co-1] = float(words[signature])
        count_co += 1
    co.close()

    p_co_copy = p_co.copy()
    i_co = 0
    while(i_co < 96):
        j_co = int(i_co/4)
        count_co = 0
        while(count_co != 4):
            for k_co in range(30):
                for l_co in range(4):
                    p_co_copy[k_co, i_co+4*count_co +
                              l_co] = p_co[k_co, j_co+l_co]
            j_co += 24
            count_co += 1
        i_co += 16
    p_co = p_co_copy.copy()

    labels = 96*[0]
    for i in range(96):
        first = i % 16
        if(first == 0 or first == 1 or first == 2 or first == 3):
            labels[i] = 'A'
        if(first == 4 or first == 5 or first == 6 or first == 7):
            labels[i] = 'C'
        if(first == 8 or first == 9 or first == 10 or first == 11):
            labels[i] = 'G'
        if(first == 12 or first == 13 or first == 14 or first == 15):
            labels[i] = 'T'
        for j in range(16):
            if(i == j):
                labels[i] += '(C>A)'
            if(i == j+16):
                labels[i] += '(C>G)'
            if(i == j+32):
                labels[i] += '(C>T)'
            if(i == j+48):
                labels[i] += '(T>A)'
            if(i == j+64):
                labels[i] += '(T>C)'
            if(i == j+80):
                labels[i] += '(T>G)'
        second = i % 4
        if(second == 0):
            labels[i] += 'A'
        if(second == 1):
            labels[i] += 'C'
        if(second == 2):
            labels[i] += 'G'
        if(second == 3):
            labels[i] += 'T'

    colorlist = 96*[0]
    for i in range(16):
        colorlist[i] = 'r'
    for i in range(16, 32):
        colorlist[i] = 'g'
    for i in range(32, 48):
        colorlist[i] = 'b'
    for i in range(48, 64):
        colorlist[i] = 'c'
    for i in range(64, 80):
        colorlist[i] = 'm'
    for i in range(80, 96):
        colorlist[i] = 'y'

    for i in range(30):
        fig = plt.figure(figsize=(6, 2))
        left = np.arange(1, 97, 1)
        height_ex = p_co[i]
        title1 = 'COSMIC Known Signature ' + str(i+1) + ' (mutation)'
        ax1 = fig.add_subplot(111)
        ax1.bar(left, height_ex, width=1, color=colorlist, align="center")
        height_limits = []
        for j in range(len(height_ex)):
            height_limits.append(height_ex[j])
        max_height_limit = max(height_limits)
        upper_lim = math.ceil(max_height_limit*10)/10
        ax1.set_ylim(0, upper_lim)
        ax1.set_xticks(left)
        ax1.set_xticklabels(labels)
        for tick in ax1.get_xticklabels():
            tick.set_rotation(90)
        ax1.tick_params(labelsize=4)
        ax1.set_xlabel('mutation x', fontsize=7)
        ax1.set_ylabel('p (mutation = x)', fontsize=7)
        ax1.set_title(title1, fontsize=7)
        fig.tight_layout()
        name = 'result/known_signatures/figure/known_' + str(i+1) + '.png'
        fig.savefig(name, dpi=200)
        plt.close(1)


if __name__ == '__main__':
    main()
