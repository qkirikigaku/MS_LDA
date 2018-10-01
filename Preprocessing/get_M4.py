import sys
import re
import pandas as pd
import numpy as np
import linecache


def main():
    args = sys.argv
    """args[1]:threshold number of word [2]:cancer type"""

    dictionary = load_dictionary()
    temp_dictionary = list()
    dictionary_indel = load_dictionary_indel()
    temp_dictionary_indel = list()

    K = 1556
    threshold = int(args[1])
    pre_data_file = 'data/data3/Pre_data3_o' + args[1] + '.txt'
    pre_data = pd.read_csv(pre_data_file, delimiter='\t', dtype={
                           'Sample name': str, 'Mutation strand': str})
    cancer_type = args[2]
    pre_data = select_cancer_type(pre_data, threshold, cancer_type)
    pre_data.reset_index(drop=True, inplace=True)

    output_file = 'data/data4_o' + args[1] + '_' + args[2] + '.txt'
    output = open(output_file, 'w')

    name_list = pre_data['Sample name']
    last_document = 'xxxxxxxxxxxx'
    number_of_document = 0
    for name in name_list:
        if(name != str(last_document)):
            last_document = name
            number_of_document += 1

    last_document = name_list[0]
    type_list = list()
    type_list.append(pre_data['Primary site'][0])
    data_mat = np.zeros([number_of_document, K])
    index = 0
    error = 0
    unexpected_error = 0
    document = 0
    for name in name_list:
        if(name != str(last_document)):
            last_document = name
            document += 1
            type_list.append(pre_data['Primary site'][index])
        CDS = list(pre_data['Mutation CDS'][index])
        if((CDS[len(CDS)-2] == '>') and
                (CDS[len(CDS)-4] not in ['A', 'C', 'G', 'T'])):
            selected = calc_word(
                pre_data['Mutation CDS'][index],
                pre_data['Mutation genome position'][index],
                pre_data['Mutation strand'][index])
            if(selected == -1):
                error += 1
            else:
                temp_dictionary = write_temp_dictionary(
                    temp_dictionary, dictionary, selected)
                data_mat[document, selected] += 1
        else:
            selected = calc_other_word(
                pre_data['Mutation CDS'][index],
                pre_data['Mutation genome position'][index],
                pre_data['Mutation strand'][index])
            if(selected not in [0, -1]):
                temp_dictionary_indel = write_temp_dictionary_indel(
                    temp_dictionary_indel, dictionary_indel, selected)
                print(temp_dictionary_indel)
                data_mat[document][selected] += 1
            elif(selected == 0):
                print('error: ' + pre_data['Mutation CDS'][index])
                error += 1
            elif(selected == -1):
                print('error: ' + pre_data['Mutation CDS'][index])
                unexpected_error += 1
        index += 1
        print(index)

    print('normal error : ' + str(error))
    print('unexpected error : ' + str(unexpected_error))

    # check
    drop_list = list()
    for i in range(number_of_document):
        sum_words = 0
        for j in range(len(temp_dictionary) + 20):
            sum_words += data_mat[i, j]
        if(sum_words == 0):
            drop_list.append(i)
    number_of_document -= len(drop_list)
    data_mat = np.delete(data_mat, drop_list, 0)
    type_list = np.delete(type_list, drop_list, 0)

    print(str(len(drop_list)))

    re_data = rewrite_data_mat(data_mat, dictionary, temp_dictionary,
                               number_of_document, dictionary_indel,
                               temp_dictionary_indel)
    number_of_vocab = len(temp_dictionary) + len(temp_dictionary_indel)
    output.write(str(number_of_document) + ' ')
    output.write(str(number_of_vocab) + '\n')

    for i in range(number_of_document):
        for j in range(number_of_vocab):
            if(j != number_of_vocab - 1):
                output.write(str(int(re_data[i][j])) + ' ')
            else:
                output.write(str(int(re_data[i][j])))
        output.write('\n')

    temp_dict_file = 'data/dictionary/data4_o' + \
        args[1] + '_' + args[2] + '.txt'
    temp_dict_indel_file = 'data/dictionary/data4_o' + \
        args[1] + '_' + args[2] + '_indel.txt'

    temp_dict = open(temp_dict_file, 'w')
    for i in range(len(temp_dictionary)):
        temp_dict.write(temp_dictionary[i] + '\n')
    temp_dict_indel = open(temp_dict_indel_file, 'w')
    for i in range(len(temp_dictionary_indel)):
        temp_dict_indel.write(temp_dictionary_indel[i] + '\n')


def select_cancer_type(pre_data, threshold, cancer_type):
    if(cancer_type != 'all'):
        pre_data = pre_data.loc[pre_data['Primary site'] == cancer_type]
    pre_data.sort_values(by='Sample name', inplace=True)
    print(pre_data)
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
    if(len(position_list) != 3):
        print('position error')
        return -1
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
    target_line = linecache.getline(GRCh_file, int(quotient)+1)

    if(((target_line[target_index] != before) and (strand == '+')) or
            ((target_line[target_index] != swap(before))and(strand == '-'))):
        print('error: ' + mutation)
        print('target: ' + target_line[target_index])
        print('strand: ' + strand)
        strand = swap(strand)
        if(((target_line[target_index] != before) and (strand == '+')) or
                ((target_line[target_index] != swap(before))and(strand == '-'))):
            print('still error')
            return -1

    if((target_index >= 2) and (target_index <= 47)):
        pattern = 1
    elif(target_index == 0):
        pattern = 2
    elif(target_index == 1):
        pattern = 3
    elif(target_index == 48):
        pattern = 4
    else:
        pattern = 5

    if(pattern == 1):
        forward = target_line[target_index - 1]
        for_forward = target_line[target_index - 2]
        backward = target_line[target_index + 1]
        back_backward = target_line[target_index + 2]
    elif(pattern == 2):
        pre_line = linecache.getline(GRCh_file, int(quotient))
        forward = pre_line[49]
        for_forward = pre_line[48]
        backward = target_line[target_index + 1]
        back_backward = target_line[target_index + 2]
    elif(pattern == 3):
        pre_line = linecache.getline(GRCh_file, int(quotient))
        for_forward = pre_line[49]
        forward = target_line[target_index - 1]
        backward = target_line[target_index + 1]
        back_backward = target_line[target_index + 2]
    elif(pattern == 4):
        post_line = linecache.getline(GRCh_file, int(quotient)+2)
        back_backward = post_line[0]
        forward = target_line[target_index - 1]
        for_forward = target_line[target_index - 2]
        backward = target_line[target_index + 1]
    if(pattern == 5):
        post_line = linecache.getline(GRCh_file, int(quotient)+2)
        backward = post_line[0]
        back_backward = post_line[1]
        forward = target_line[target_index - 1]
        for_forward = target_line[target_index - 2]

    if(((strand == '+') and (before in ['A', 'G'])) or
            ((strand == '-') and (before in ['C', 'T']))):
        buf_f = swap(forward)
        buf_ff = swap(for_forward)
        forward = swap(backward)
        for_forward = swap(back_backward)
        backward = buf_f
        back_backward = buf_ff
    if(before in ['A', 'G']):
        before = swap(before)
        after = swap(after)

    if(for_forward == 'A'):
        first = 0
    elif(for_forward == 'C'):
        first = 1
    elif(for_forward == 'G'):
        first = 2
    else:
        first = 3

    if(forward == 'A'):
        second = 0
    elif(forward == 'C'):
        second = 1
    elif(forward == 'G'):
        second = 2
    else:
        second = 3

    if(before == 'C'):
        if(after == 'A'):
            third = 0
        elif(after == 'G'):
            third = 1
        else:
            third = 2
    elif(before == 'T'):
        if(after == 'A'):
            third = 3
        elif(after == 'C'):
            third = 4
        else:
            third = 5
    elif(before == 'G'):
        if(after == 'T'):
            third = 0
        elif(after == 'C'):
            third = 1
        else:
            third = 2
    else:
        if(after == 'T'):
            third = 3
        elif(after == 'G'):
            third = 4
        else:
            third = 5

    if(back_backward == 'A'):
        fourth = 0
    elif(back_backward == 'C'):
        fourth = 1
    elif(back_backward == 'G'):
        fourth = 2
    else:
        fourth = 3

    if(backward == 'A'):
        fifth = 0
    elif(backward == 'C'):
        fifth = 1
    elif(backward == 'G'):
        fifth = 2
    else:
        fifth = 3
    answer = 384*first + 96*second + 16*third + 4*fourth + fifth
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


def load_dictionary_indel():
    dictionary_file = 'data/dictionary_indel.txt'
    file = open(dictionary_file, 'r')
    lines = file.readlines()
    dictionary = list()
    for line in lines:
        text = line[:-1]
        dictionary.append(text)
    return dictionary


def write_temp_dictionary(temp_dictionary, dictionary, pre_index):
    flag = 0
    for i in temp_dictionary:
        if(i == dictionary[pre_index]):
            flag = 1
            break
    if(flag == 0):
        temp_dictionary.append(dictionary[pre_index])
    return temp_dictionary


def write_temp_dictionary_indel(temp_dictionary_indel,
                                dictionary_indel, pre_index):
    flag = 0
    index = pre_index - 1536
    for i in temp_dictionary_indel:
        if(i == dictionary_indel[index]):
            flag = 1
            break
    if(flag == 0):
        temp_dictionary_indel.append(dictionary_indel[index])
    return temp_dictionary_indel


def rewrite_data_mat(data_mat, dictionary, temp_dictionary,
                     number_of_document, dictionary_indel, temp_dictionary_indel):
    number_of_vocab = len(temp_dictionary) + len(temp_dictionary_indel)
    re_data = [[0 for i in range(number_of_vocab)]
               for j in range(number_of_document)]
    for i in range(number_of_document):
        for j in range(1556):
            if((data_mat[i, j] != 0) and (j < 1536)):
                index = temp_dictionary.index(dictionary[j])
                re_data[i][index] = data_mat[i, j]
            elif((data_mat[i, j] != 0) and (j >= 1536)):
                index = temp_dictionary_indel.index(
                    dictionary_indel[j-1536]) + len(temp_dictionary)
                re_data[i][index] = data_mat[i, j]
    return re_data


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
            return 1536
        else:
            return 1546
    elif(new_vocab == 'CG'):
        if(num >= 10):
            return 1537
        else:
            return 1547
    elif(new_vocab == 'GC'):
        if(num >= 10):
            return 1538
        else:
            return 1548
    elif(new_vocab == 'TA'):
        if(num >= 10):
            return 1539
        else:
            return 1549
    elif(new_vocab == 'AC'):
        if(num >= 10):
            return 1540
        else:
            return 1550
    elif(new_vocab == 'CT'):
        if(num >= 10):
            return 1541
        else:
            return 1551
    elif(new_vocab == 'TG'):
        if(num >= 10):
            return 1542
        else:
            return 1552
    elif(new_vocab == 'GA'):
        if(num >= 10):
            return 1543
        else:
            return 1553
    elif(new_vocab == 'AA'):
        if(num >= 10):
            return 1544
        else:
            return 1554
    elif(new_vocab == 'CC'):
        if(num >= 10):
            return 1545
        else:
            return 1555
    else:
        print(new_vocab + '\n')
        return -1


if __name__ == '__main__':
    main()
