import sys
import numpy as np
import pandas as pd

def main():
    args = sys.argv
    """ args[1]:data_type [2]:threshold number word """
    threshold = int(args[2])
    
    input_file = 'data/data' + args[1] + '/Pre_data' + args[1] + '_o' + args[2] + '_hist.txt'
    file = pd.read_csv(input_file, delimiter = '\t')
    data = file.copy()
    del file
    check_whole_document(data)
    print('------------------------------------------------------------')
    check_Primary_Histology(data)
    print('------------------------------------------------------------')
    check_Histology_subtype(data)

def check_whole_document(data):
    temp = data.copy()
    document_list = list(temp['Sample name'])
    number = 1
    last_document = document_list[0]
    for document in document_list:
        if(document != last_document):
            last_document = document
            number += 1
    print('the number of document is ...')
    print(str(number-1) + '\n')

def check_Primary_Histology(data):
    temp = data.copy()
    del data
    temp.sort_values(by=['Primary histology','Sample name'], inplace=True)
    temp.reset_index(drop=True, inplace=True)
    last_name  = temp['Sample name'][0] 
    drop_list = list()
    name_list = list(temp['Sample name'])
    index = 0
    for name in name_list:
        if(name != last_name):
            last_name = name
        else:
            drop_list.append(index)
        index += 1
    temp.drop(drop_list, inplace=True)
    
    List = list()
    site_number_list = [0]
    site_list = list(temp['Primary histology'])
    hist_list = list(temp['Primary site'])
    hist_type_list = [[]]
    List.append(site_list[0])
    last_site = site_list[0]
    index = 0
    for i, site in enumerate(site_list):
        if(site != last_site):
            last_site = site
            List.append(site)
            site_number_list.append(0)
            hist_type_list.append([])
            index += 1
        site_number_list[index] += 1
        hist_type_list[index].append(hist_list[i])
    
    site_num = len(site_number_list)
    hist_dict_list = []
    for i in range(site_num):
        hist_dict_list.append({})
        for x in hist_type_list[i]:
            hist_dict_list[i][x] = hist_type_list[i].count(x)

    for i, site in enumerate(List):
        string = site + ' :'
        print(string + str(site_number_list[i]))
        for key,val in sorted(hist_dict_list[i].items(), key=lambda x: -x[1]):
            print('∟ ' + str(key) + ' : ' + str(val))
        print('\n')

def check_Histology_subtype(data):
    temp = data.copy()
    del data
    temp.sort_values(by=['Histology subtype 1','Sample name'], inplace=True)
    temp.reset_index(drop=True, inplace=True)
    last_name  = temp['Sample name'][0] 
    drop_list = list()
    name_list = list(temp['Sample name'])
    index = 0
    for name in name_list:
        if(name != last_name):
            last_name = name
        else:
            drop_list.append(index)
        index += 1
    temp.drop(drop_list, inplace=True)
    
    List = list()
    site_number_list = [0]
    site_list = list(temp['Histology subtype 1'])
    hist_list = list(temp['Primary site'])
    hist_type_list = [[]]
    List.append(site_list[0])
    last_site = site_list[0]
    index = 0
    for i, site in enumerate(site_list):
        if(site != last_site):
            last_site = site
            List.append(site)
            site_number_list.append(0)
            hist_type_list.append([])
            index += 1
        site_number_list[index] += 1
        hist_type_list[index].append(hist_list[i])
    
    site_num = len(site_number_list)
    hist_dict_list = []
    for i in range(site_num):
        hist_dict_list.append({})
        for x in hist_type_list[i]:
            hist_dict_list[i][x] = hist_type_list[i].count(x)

    index = 0
    for i, site in enumerate(List):
        string = site + ' :'
        print(string + str(site_number_list[i]))
        for key,val in sorted(hist_dict_list[i].items(), key=lambda x: -x[1]):
            print('∟ ' + str(key) + ' : ' + str(val))
        print('\n')


main()
