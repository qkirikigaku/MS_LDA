import pandas as pd
import re
import linecache


def main():
    mutation_file = 'raw_data/CosmicMutantExport.tsv'
    complete_data = pd.read_csv(mutation_file, delimiter='\t')
    make_dictionary(complete_data)


def make_dictionary(complete_data):
    dictionary_file = 'data/dictionary.txt'
    dictionary = open(dictionary_file, 'w')
    dictionary = extract_base_substitution(dictionary)
    dictionary = extract_multi_substitution(dictionary, complete_data)
    dictionary = extract_indel(dictionary, complete_data)


def extract_base_substitution(dictionary):
    # extract base substitution
    for i in range(1536):
        first = i // 384
        second = (i % 384) // 96
        third = ((i % 384) % 96) // 16
        fourth = (((i % 384) % 96) % 16) // 4
        fifth = (((i % 384) % 96) % 16) % 4
        if(first == 0):
            word = 'A'
        elif(first == 1):
            word = 'C'
        elif(first == 2):
            word = 'G'
        else:
            word = 'T'
        if(second == 0):
            word += 'A'
        elif(second == 1):
            word += 'C'
        elif(second == 2):
            word += 'G'
        else:
            word += 'T'
        if(third == 0):
            word += '[C>A]'
        elif(third == 1):
            word += '[C>G]'
        elif(third == 2):
            word += '[C>T]'
        elif(third == 3):
            word += '[T>A]'
        elif(third == 4):
            word += '[T>C]'
        else:
            word += '[T>G]'
        if(fourth == 0):
            word += 'A'
        elif(fourth == 1):
            word += 'C'
        elif(fourth == 2):
            word += 'G'
        else:
            word += 'T'
        if(fifth == 0):
            word += 'A'
        elif(fifth == 1):
            word += 'C'
        elif(fifth == 2):
            word += 'G'
        else:
            word += 'T'
        dictionary.write(word + '\n')
    print('complete base substitution...')
    return dictionary


def extract_multi_substitution(dictionary, complete_data):
    # extract multiple substitution
    data = complete_data[complete_data['Mutation Description'].
                         str.contains('Substitution')]
    data = data.sort_values(by='Mutation CDS')
    data.reset_index(drop=True, inplace=True)

    CDS_list = data['Mutation CDS']
    drop_list = list()
    index = 0
    for CDS in CDS_list:
        char_list = list(CDS)
        if((char_list[len(char_list)-4] not in ['A', 'C', 'G', 'T', '>'])
                or char_list[len(char_list)-1] == '?'):
            drop_list.append(index)
        index += 1
    data.drop(drop_list, inplace=True)
    CDS_list = data['Mutation CDS']
    multi_list = list()
    last_CDS = 'xxxxxxx'
    for CDS in CDS_list:
        if(CDS != last_CDS):
            char_list = list(CDS)
            number = 0
            for i in char_list:
                if(i in ['A', 'C', 'G', 'T']):
                    break
                number += 1
            M_index = 0
            for i in char_list:
                if(i == '>'):
                    break
                M_index += 1
            mutation_number = len(char_list) - number
            words = list()
            for i in range(mutation_number):
                words.append(char_list[i + number])
            out = ''
            flag = 0
            for word in words:
                if(word not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']):
                    out += word
                if((word == 'M')or(word == 'R')or(word == 'Y')or(word == 'K')
                        or(word == 'H')or(word == 'B')):
                    flag = 1
            out_re = ''
            for i in range(M_index - number):
                word = words[M_index - number - 1 - i]
                if(word not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']):
                    word_re = swap(word)
                    out_re += word_re
            out_re += '>'
            for i in range(len(char_list) - M_index - 1):
                word = words[len(words) - i - 1]
                if(word not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']):
                    word_re = swap(word)
                    out_re += word_re
            if(out == ''):
                flag = 1
            for i in multi_list:
                if((out == i) or (out_re == i)):
                    flag = 1
                    break
            if(flag == 0):
                multi_list.append(out)
                dictionary.write(out + '\n')
        last_CDS = CDS
    print('complete multi substitution...')
    return dictionary


def extract_indel(dictionary, complete_data):
    # extract insertion
    data = complete_data[complete_data['Mutation Description']
                         .str.contains('Insertion')]
    data = data.sort_values(by='Mutation CDS')
    data.reset_index(drop=True, inplace=True)

    CDS_list = data['Mutation CDS']
    last_CDS = 'xxxxxxxx'
    ins_list = list()
    index = 0
    for CDS in CDS_list:
        if(CDS != last_CDS):
            char_list = list(CDS)
            number = 0
            for char in char_list:
                if(char == 's'):
                    break
                number += 1
            words = list()
            for i in range(len(char_list) - number - 1):
                words.append(char_list[i + number + 1])
            out = ''
            flag = 0
            for word in words:
                if(word not in ['?', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']):
                    out += word
                if(word == 'N'):
                    flag = 1
            out_re = ''
            for i in range(len(words)):
                word = words[len(words) - i - 1]
                if(word not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']):
                    word_re = swap(word)
                    out_re += word_re
            if(out == ''):
                flag = 1
            mat = get_base('error', data['Mutation genome position'][index],
                           data['Mutation strand'][index], 'ins')
            if(mat[0] == 'error'):
                flag = 1
            five = mat[0]
            three = mat[1]
            five_re = swap(three)
            three_re = swap(five)
            out = five + '[' + out + ']' + three
            out_re = five_re + '[' + out_re + ']' + three_re
            for i in ins_list:
                if((out == i) or (out_re == i)):
                    flag = 1
                    break
            if(flag == 0):
                ins_list.append(out)
                dictionary.write('ins' + out + '\n')
            last_CDS = CDS
        index += 1
    print('complete insertion...')

    # extract deletion
    data = complete_data[complete_data['Mutation Description']
                         .str.contains('Deletion')]
    data = data.sort_values(by='Mutation CDS')
    data.reset_index(drop=True, inplace=True)

    CDS_list = data['Mutation CDS']
    last_CDS = 'xxxxxxxx'
    del_list = list()
    index = 0
    for CDS in CDS_list:
        if(CDS != last_CDS):
            char_list = list(CDS)
            number = 0
            for char in char_list:
                if(char == 'l'):
                    break
                number += 1
            words = list()
            for i in range(len(char_list) - number - 1):
                words.append(char_list[i + number + 1])
            out = ''
            flag = 0
            for word in words:
                if(word not in ['?', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']):
                    out += word
                if(word == 'N'):
                    flag = 1
            out_re = ''
            for i in range(len(words)):
                word = words[len(words) - i - 1]
                if(word not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']):
                    word_re = swap(word)
                    out_re += word_re
            if(out == ''):
                flag = 1
                top = 'error'
            else:
                top = out[0]
            mat = get_base(top, data['Mutation genome position'][index],
                           data['Mutation strand'][index], 'del')
            if(mat[0] == 'error'):
                flag = 1
            five = mat[0]
            three = mat[1]
            five_re = swap(three)
            three_re = swap(five)
            out = five + '[' + out + ']' + three
            out_re = five_re + '[' + out_re + ']' + three_re
            for i in del_list:
                if((out == i) or (out_re == i)):
                    flag = 1
                    break
            if(flag == 0):
                del_list.append(out)
                dictionary.write('del' + out + '\n')
            last_CDS = CDS
        index += 1
    print('complete deletion...')


def swap(word):
    if(word == 'A'):
        word_re = 'T'
    elif(word == 'C'):
        word_re = 'G'
    elif(word == 'G'):
        word_re = 'C'
    elif(word == 'T'):
        word_re = 'A'
    else:
        word_re = word
    return word_re


def get_base(top, position, strand, indel):
    mat = [0 for i in range(2)]
    if(position != position):
        mat[0] = 'error'
        mat[1] = 'error'
        return mat
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
    if(indel == 'ins'):
        num = 0
    if(num >= 50):
        mat[0] = 'error'
        mat[1] = 'error'
        return mat
    GRCh_file = 'raw_data/chr' + str(chromosome) + '.fa'
    quotient = start // 50
    surplus = start % 50
    if(surplus != 0):
        if(indel == 'del'):
            target_index = int(surplus) - 1
        elif(indel == 'ins'):
            target_index = int(surplus)
    elif(surplus == 0):
        if(indel == 'del'):
            quotient -= 1
            target_index = 49
        elif(indel == 'ins'):
            target_index = 0
    if((target_index > 0) and (target_index + num < 50)):
        pattern = 1
    elif(target_index == 0):
        pattern = 2
    else:
        pattern = 3
    target_line = linecache.getline(GRCh_file, int(quotient)+1)

    if(top != 'error'):
        target_base = target_line[target_index]
        if((strand == '+')and(target_base != top)):
            print('target_base:' + target_base + '  top: ' + top
                  + '  strand: ' + strand + '  target_index:'
                    + str(target_index))
            print('error!')
            mat[0] = 'error'
            mat[1] = 'error'
            return mat
        if(strand == '-'):
            if(target_index + num - 1 < 50):
                target_base = target_line[target_index + num - 1]
                if(swap(target_base) != top):
                    print('target_base:' + target_base + '  top: ' + top
                          + '  strand: ' + strand + '  target_index:'
                            + str(target_index) + '  num: ' + str(num))
                    print('error!')
                    mat[0] = 'error'
                    mat[1] = 'error'
                    return mat
            else:
                post_line = linecache.getline(GRCh_file, int(quotient)+2)
                if(swap(post_line[target_index + num - 51]) != top):
                    print('target_base:' + target_base + '  top: ' + top
                          + '  strand: ' + strand + '  target_index:'
                            + str(target_index) + '  num: ' + str(num))
                    print('error!')
                    mat[0] = 'error'
                    mat[1] = 'error'
                    return mat

    if(pattern == 1):
        mat[0] = target_line[target_index - 1]
        mat[1] = target_line[target_index + num]
    elif(pattern == 2):
        pre_line = linecache.getline(GRCh_file, int(quotient))
        mat[0] = pre_line[49]
        mat[1] = target_line[target_index + num]
    else:
        post_line = linecache.getline(GRCh_file, int(quotient)+2)
        mat[0] = target_line[target_index - 1]
        mat[1] = post_line[target_index + num - 50]
    if(strand == '-'):
        buf = swap(mat[0])
        mat[0] = swap(mat[1])
        mat[1] = buf
    return mat


if __name__ == '__main__':
    main()
