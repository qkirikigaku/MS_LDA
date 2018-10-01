import sys
import re
import pandas as pd
import numpy as np
import linecache


def main():
    args = sys.argv
    """args[1] : threshold number of word [2]:cancer type"""

    dictionary_indel = load_dictionary_indel()
    temp_dictionary_indel = list()

    K = 116
    threshold = int(args[1])
    pre_data_file = 'data/data3/Pre_data3_o' + args[1] + '.txt'
    pre_data = pd.read_csv(pre_data_file, delimiter='\t', dtype={
                           'Sample name': str, 'Mutation strand': str})
    cancer_type = args[2]
    pre_data = select_cancer_type(pre_data, threshold, cancer_type)
    print(pre_data)
    pre_data.reset_index(drop=True, inplace=True)

    output_file = 'data3_o' + args[1] + '_' + args[2] + '.txt'
    output = open('data/' + output_file, 'w')

    name_list = pre_data['Sample name']
    last_document = 'xxxxxxxxxxxx'
    number_of_document = 0
    for name in name_list:
        if(name != last_document):
            last_document = name
            number_of_document += 1

    data_mat = [[0 for i in range(K)] for j in range(number_of_document)]

    document = 0
    index = 0
    last_document = name_list[0]
    error = 0
    unexpected_error = 0
    for name in name_list:
        if(name != last_document):
            last_document = name
            document += 1
        CDS = list(pre_data['Mutation CDS'][index])
        if((CDS[len(CDS)-2] == '>') and 
                (CDS[len(CDS)-4] not in ['A', 'C', 'G', 'T'])):
            selected = calc_word(
                pre_data['Mutation CDS'][index], 
                pre_data['Mutation genome position'][index], 
                pre_data['Mutation strand'][index])
            if(selected != -1):
                data_mat[document][selected] += 1
        else:
            selected = calc_other_word(
                pre_data['Mutation CDS'][index], 
                pre_data['Mutation genome position'][index], 
                pre_data['Mutation strand'][index])
            if(selected not in [0, -1]):
                temp_dictionary_indel = write_temp_dictionary_indel(
                    temp_dictionary_indel, dictionary_indel, selected)
                data_mat[document][selected] += 1
            elif(selected == 0):
                print('error: ' + pre_data['Mutation CDS'][index])
                error += 1
            elif(selected == -1):
                print('error: ' + pre_data['Mutation CDS'][index])
                unexpected_error += 1
        index += 1
    print('normal error : ' + str(error))
    print('unexpected error : ' + str(unexpected_error))

    re_data = rewrite_data_mat(
        data_mat, number_of_document, dictionary_indel, temp_dictionary_indel)
    number_of_vocab = 96 + len(temp_dictionary_indel)

    output.write(str(number_of_document) + ' ')
    output.write(str(number_of_vocab) + '\n')
    for i in range(number_of_document):
        for j in range(number_of_vocab):
            if(j != number_of_vocab - 1):
                output.write(str(re_data[i][j]) + ' ')
            else:
                output.write(str(re_data[i][j]))
        output.write('\n')

    temp_dict_indel_file = 'data/dictionary/data3_o' + \
        args[1] + '_' + args[2] + '_indel.txt'
    temp_dict_indel = open(temp_dict_indel_file, 'w')
    for i in range(len(temp_dictionary_indel)):
        temp_dict_indel.write(temp_dictionary_indel[i] + '\n')


def select_cancer_type(pre_data, threshold, cancer_type):
    pre_data = pre_data.loc[pre_data['Primary site'] == cancer_type]
    pre_data.sort_values(by='Sample name', inplace=True)
    pre_data.reset_index(drop=True, inplace=True)
    name_list = pre_data['Sample name']
    last_document = name_list[0]
    drop_list = list()
    sum_of_words = 0
    temp_index = 0
    index_list = list()
    for name in name_list:
        if(name != last_document):
            if(sum_of_words < threshold):
                drop_list.extend(index_list)
            sum_of_words = 0
            last_document = name
            index_list = list()
        sum_of_words += 1
        index_list.append(temp_index)
        temp_index += 1
    if(sum_of_words < threshold):
        drop_list.extend(index_list)
    pre_data.drop(drop_list, inplace=True)
    return pre_data


def calc_word(mutation, position, strand):
    before = mutation[len(mutation)-3]
    after = mutation[len(mutation)-1]
    position_list = re.split(r'[:-]', position)
    if(int(position_list[0]) == 23):
        chromosome = 'X'
    elif(int(position_list[0]) == 24):
        chromosome = 'Y'
    elif(int(position_list[0]) == 25):
        chromosome = 'M'
    else:
        chromosome = int(position_list[0])
    start = int(position_list[1])
    num = int(position_list[2]) - int(position_list[1]) + 1
    GRCh_file = 'raw_data/chr' + str(chromosome) + '.fa'
    quotient = start // 50
    surplus = start % 50

    if(surplus != 0):
        target_index = int(surplus) - 1
    else:
        quotient -= 1
        target_index = 49
    targetline = linecache.getline(GRCh_file, int(quotient)+1)

    if(((targetline[target_index] != before) and (strand == '+')) or 
            ((targetline[target_index] != swap(before))and(strand == '-'))):
        print('error: ' + mutation)
        print('target: ' + targetline[target_index])
        print('strand: ' + strand)
        strand = swap(strand)
        if(((targetline[target_index] != before) and (strand == '+')) or 
                ((targetline[target_index] != swap(before))and(strand == '-'))):
            print('still error')
            return -1

    if((target_index >= 1) and (target_index <= 48)):
        pattern = 1
    elif(target_index == 0):
        pattern = 2
    elif(target_index == 49):
        pattern = 3

    if(pattern == 1):
        forward = targetline[target_index - 1]
        backward = targetline[target_index + 1]
    elif(pattern == 2):
        pre_line = linecache.getline(GRCh_file, int(quotient))
        forward = pre_line[49]
        backward = targetline[target_index + 1]
    elif(pattern == 3):
        post_line = linecache.getline(GRCh_file, int(quotient)+2)
        forward = targetline[target_index - 1]
        backward = post_line[0]

    if(((strand == '+') and (before in ['A', 'G'])) or 
            ((strand == '-') and (before in ['C', 'T']))):
        buf_f = swap(forward)
        forward = swap(backward)
        backward = buf_f
    if(before in ['A', 'G']):
        before = swap(before)
        after = swap(after)

    if(forward == 'A'):
        first = 0
    elif(forward == 'C'):
        first = 1
    elif(forward == 'G'):
        first = 2
    else:
        first = 3

    if(before == 'C'):
        if(after == 'A'):
            second = 0
        elif(after == 'G'):
            second = 1
        else:
            second = 2
    elif(before == 'T'):
        if(after == 'A'):
            second = 3
        elif(after == 'C'):
            second = 4
        else:
            second = 5

    if(backward == 'A'):
        third = 0
    elif(backward == 'C'):
        third = 1
    elif(backward == 'G'):
        third = 2
    else:
        third = 3
    answer = 24*first + 4*second + third
    return(answer)


def swap(base):
    if(base == 'A'):
        return('T')
    elif(base == 'C'):
        return('G')
    elif(base == 'G'):
        return('C')
    elif(base == 'T'):
        return('A')
    elif(base == '+'):
        return('-')
    elif(base == '-'):
        return('+')
    else:
        return(base)


def load_dictionary():
    dictionary_file = 'data/dictionary.txt'
    file = open(dictionary_file, 'r')
    lines = file.readlines()
    dictionary = list()
    for line in lines:
        text = line[:-1]
        dictionary.append(text)
    return dictionary


def calc_other_word(mutation, position, strand):
    if(mutation.find('>') > -1):
        return 0
    else:
        words = list(mutation)
        if(mutation.find('ins') > -1):
            number = mutation.find('s')
            head = 'ins'
        elif(mutation.find('del') > -1):
            number = mutation.find('l')
            head = 'del'
        else:
            print('errorだよ')
            return 0

    if(words[number + 1] in ['1', '2', '3', '4', '5', '6', '7', '8', '9']):
        leng = ''
        for i in range(number + 1, len(words)):
            leng += words[i]
        length = int(leng)
    elif(words[number + 1] in ['A', 'C', 'G', 'T']):
        vocab = list()
        for i in range(number + 1, len(words)):
            vocab += words[i]
        length = len(vocab)
    else:
        return 0

    position_list = re.split(r'[:-]', position)
    if(int(position_list[0]) == 23):
        chromosome = 'X'
    elif(int(position_list[0]) == 24):
        chromosome = 'Y'
    elif(int(position_list[0]) == 25):
        chromosome = 'M'
    else:
        chromosome = int(position_list[0])
    start = int(position_list[1])
    end = int(position_list[2])
    if(head == 'ins'):
        num = 0
    elif(head == 'del'):
        num = length
    GRCh_file = 'raw_data/chr' + str(chromosome) + '.fa'
    quotient = start // 50
    surplus = start % 50
    if(head == 'ins'):
        if(surplus != 0):
            forward_index = int(surplus) - 1
        else:
            quotient -= 1
            forward_index = 49
    else:
        if(surplus not in [0, 1]):
            forward_index = int(surplus) - 2
        else:
            quotient -= 1
            forward_index = int(surplus) + 48
    targetline = linecache.getline(GRCh_file, int(quotient)+1)
    forward = targetline[forward_index]
    if(head == 'ins'):
        if(forward_index != 49):
            backward_index = forward_index + 1
        else:
            quotient += 1
            backward_index = 0
        targetline = linecache.getline(GRCh_file, int(quotient) + 1)
        backward = targetline[backward_index]
    else:
        end_quotient = end // 50
        end_surplus = end % 50
        if(end_surplus != 0):
            backward_index = end_surplus
        else:
            backward_index = end_surplus
            end_quotient += 1
        endline = linecache.getline(GRCh_file, int(end_quotient) + 1)
        backward = endline[backward_index]
    vocab = forward + backward
    return make_vocab(vocab, length)


def make_vocab(vocab, num):
    if(vocab in ['GT', 'AG', 'CA', 'TC', 'TT', 'GG']):
        new_vocab = swap(vocab[1]) + swap(vocab[0])
    else:
        new_vocab = vocab
    if(new_vocab == 'AT'):
        if(num >= 10):
            return 96
        else:
            return 106
    elif(new_vocab == 'CG'):
        if(num >= 10):
            return 97
        else:
            return 107
    elif(new_vocab == 'GC'):
        if(num >= 10):
            return 98
        else:
            return 108
    elif(new_vocab == 'TA'):
        if(num >= 10):
            return 99
        else:
            return 109
    elif(new_vocab == 'AC'):
        if(num >= 10):
            return 100
        else:
            return 110
    elif(new_vocab == 'CT'):
        if(num >= 10):
            return 101
        else:
            return 111
    elif(new_vocab == 'TG'):
        if(num >= 10):
            return 102
        else:
            return 112
    elif(new_vocab == 'GA'):
        if(num >= 10):
            return 103
        else:
            return 113
    elif(new_vocab == 'AA'):
        if(num >= 10):
            return 104
        else:
            return 114
    elif(new_vocab == 'CC'):
        if(num >= 10):
            return 105
        else:
            return 115
    else:
        print(new_vocab + '\n')
        return -1


def load_dictionary_indel():
    dictionary_file = 'data/dictionary_indel.txt'
    file = open(dictionary_file, 'r')
    lines = file.readlines()
    dictionary = list()
    for line in lines:
        text = line[:-1]
        dictionary.append(text)
    return dictionary


def write_temp_dictionary_indel(temp_dictionary_indel, 
        dictionary_indel, pre_index):
    flag = 0
    index = pre_index - 96
    for i in temp_dictionary_indel:
        if(i == dictionary_indel[index]):
            flag = 1
            break
    if(flag == 0):
        temp_dictionary_indel.append(dictionary_indel[index])
    return temp_dictionary_indel


def rewrite_data_mat(data_mat, number_of_document, 
        dictionary_indel, temp_dictionary_indel):
    number_of_vocab = 96 + len(temp_dictionary_indel)
    re_data = [[0 for i in range(number_of_vocab)]
               for j in range(number_of_document)]
    for i in range(number_of_document):
        for j in range(116):
            if(j < 96):
                re_data[i][j] = data_mat[i][j]
            elif((data_mat[i][j] != 0) and (j >= 96)):
                index = temp_dictionary_indel.index(
                    dictionary_indel[j-96]) + 96
                re_data[i][index] = data_mat[i][j]
    return re_data


if __name__ == '__main__':
    main()
