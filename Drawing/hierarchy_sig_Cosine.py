import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn import manifold
from mpl_toolkits.mplot3d.axes3d import Axes3D
import sys
import n2_to_n1 as n2
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram

def main():
    args = sys.argv
    #[1] : data_type [2] : threshold typeは30文書以上あるデータを使用
    
    types = ['breast', 'endometrium', 'large_intestine', 'liver', 'lung'] \
            + ['oesophagus', 'prostate', 'skin', 'soft_tissue', 'stomach'] \
            + ['upper_aerodigestive_tract', 'urinary_tract']
    
    data_type = args[1]
    if(data_type == '1'):
        V = 96
        num_topics = [7, 12, 18, 11, 6, 7, 8, 27, 4, 6, 13, 5]
    elif(data_type == '2'):
        V = 1536
        num_topics = [6, 7, 13, 10, 5, 6, 4, 10, 4, 7, 5, 5]
    elif(data_type == '3'):
        V = 116
        num_topics = [10, 13, 22, 13, 7, 6, 7, 25, 5, 14, 12, 5]
    elif(data_type == '4'):
        V = 1556
        num_topics = [5, 8, 12, 9, 6, 6, 4, 13, 4, 11, 5, 5]

    threshold = args[2]
    
    topic_mat = load_result(data_type, threshold, types, num_topics, V)
    
    if(data_type == '1'):
        types.append('Known')
        num_topics.append(30)
        topic_mat = add_Known(topic_mat)


    do_hierarchy(data_type, threshold, types, num_topics, V, topic_mat)

def load_result(data_type, threshold, types, num_topics, V):
    num_all_topic = 0
    for i in range(len(num_topics)):
        num_all_topic += num_topics[i]
    topic_mat = np.zeros([num_all_topic, V])
    
    temp_sig = 0
    for i in range(len(types)):
        file_name = 'result/data' + data_type + '_o' + threshold + '_' \
                + types[i] + '/result_k'
        if(num_topics[i] < 10):
            file_name += '0' + str(num_topics[i]) + '.txt'
        else:
            file_name += str(num_topics[i]) + '.txt'
        lines = open(file_name).readlines()

        if(data_type == '1'):
            temp_dictionary = []
            dictionary = []
            temp_dictionary_indel = []
            dictionary_indel = []
        elif(data_type == '2'):
            dictionary_file = 'data/dictionary/data2_o' + threshold + '_' \
                    + types[i] + '.txt'
            temp_dictionary = n2.load_dictionary(dictionary_file)
            dictionary = n2.load_dictionary('data/dictionary.txt')
            temp_dictionary_indel = []
            dictionary_indel = []
        elif(data_type == '3'):
            temp_dictionary = []
            dictionary = []
            dictionary_file_indel = 'data/dictionary/data3_o' + threshold \
                    + '_' + types[i] + '_indel.txt'
            temp_dictionary_indel = n2.load_dictionary(dictionary_file_indel)
            dictionary_indel = n2.load_dictionary('data/dictionary_indel.txt')
        elif(data_type == '4'):
            dictionary_file = 'data/dictionary/data4_o' + threshold \
                    + '_' + types[i] + '.txt'
            temp_dictionary = n2.load_dictionary(dictionary_file)
            dictionary = n2.load_dictionary('data/dictionary.txt')
            dictionary_file_indel = 'data/dictionary/data4_o' + threshold \
                    + '_' + types[i] + '_indel.txt'
            temp_dictionary_indel = n2.load_dictionary(dictionary_file_indel)
            dictionary_indel = n2.load_dictionary('data/dictionary_indel.txt')
        

        count = 0
        for line in lines:
            if((count > 1) and (count <= 1 + num_topics[i])):
                signatures = line.split()
                for v in range(len(signatures)):
                    index = ref_dic(data_type, v, temp_dictionary, dictionary, \
                            temp_dictionary_indel, dictionary_indel)
                    topic_mat[temp_sig, index] = float(signatures[v])
                temp_sig += 1
            count += 1
            
    return topic_mat

def add_Known(topic_mat):
    new_mat = np.zeros([len(topic_mat) + 30, 96])
    for i in range(len(topic_mat)):
        for j in range(96):
            new_mat[i,j] = topic_mat[i,j]
    lines = open('data/signature_probability.txt').readlines()
    V = 96
    count = 0
    temp_sig = 0
    for line in lines:
        if(count != 0):
            signatures = line.split()
            for i in range(30):
                new_mat[len(new_mat) - 30 + i, count - 1] \
                        = float(signatures[i + 3])
            temp_sig += 1
        count += 1
    return new_mat


def ref_dic(data_type, v, temp_dictionary, dictionary, \
        temp_dictionary_indel, dictionary_indel):
    if(data_type == '1'):
        return v
    elif(data_type == '2'):
        return dictionary.index(temp_dictionary[v])
    elif(data_type == '3'):
        if(v < 96):
            return v
        else:
            return 96 + dictionary_indel.index(temp_dictionary_indel[v-96])
    elif(data_type == '4'):
        if(v < len(temp_dictionary)):
            return dictionary.index(temp_dictionary[v])
        else:
            return 1536 + dictionary_indel.index(temp_dictionary_indel[v \
                    - len(temp_dictionary)])

def do_hierarchy(data_type, threshold, types, num_topics, V, topic_mat):
    dist_mat = pdist(topic_mat, 'cosine')
    result = linkage(dist_mat, method = 'average')
    

    name = []
    index = 0
    for i in range(len(types)):
        for j in range(num_topics[i]):
            name.append(str(j+1) + ' ')
            name[index] += '(' + types[i] + ')'
            index += 1

    fig = plt.figure(figsize=(12,12))
    result = dendrogram(result, labels = name,
               orientation = 'right',
               color_threshold = 0.2,
               above_threshold_color = '#BCBDDC')

    print(len(result['icoord']))
    print(len(result['dcoord']))
    print(len(result['ivl']))
    print(result['color_list'])

    file_name = 'result/general/data' + data_type + '_o' \
            + threshold + '_hierarchy_Cosine.png'
    fig.tight_layout()
    
    plt.savefig(file_name, dpi=300)

def select_type(i, num_topics):
    x = i
    for j in range(len(num_topics)):
        x -= num_topics[j]
        if(x < 0):
            return j

def select_types(i, num_topics):
    x = i
    for j in range(len(num_topics)):
        x -= num_topics[j]
        if(x < 0):
            x += num_topics[j]
            return x+1


if __name__ == '__main__':
    main()
