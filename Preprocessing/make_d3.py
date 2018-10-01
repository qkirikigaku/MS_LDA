import sys
import pandas as pd


def main():
    args = sys.argv
    """ args[1]:array_number [2]:division number"""
    index = int(args[1])
    Num = int(args[2])
    mutation_file = 'raw_data/CosmicMutantExport.tsv'
    complete_data = pd.read_csv(mutation_file, delimiter='\t')
    data = select_data(complete_data, index, Num)
    output_file = 'data/data3/Pre_data3_' + args[1] + '.txt'
    data.to_csv(output_file, sep='\t', index=None)


def select_data(complete_data, index, Num):
    """ make_data 3 """
    data = complete_data[complete_data['Mutation Description']\
                         .str.contains('Substitution')\
                         | complete_data['Mutation Description']\
                         .str.contains('Insertion')\
                         | complete_data['Mutation Description']\
                         .str.contains('Deletion')]
    del complete_data

    """ select which use GRCh38 """
    data = data[data['GRCh'] == 38]

    data = data.loc[:, ['Sample name', 'Mutation CDS', 'Mutation Description',
                        'Mutation genome position', 'Primary site',
                        'Primary histology', 'Histology subtype 1',
                        'Mutation strand']]
    data = data.sort_values(by='Sample name')
    data.reset_index(drop=True, inplace=True)

    """ Pearallelization """
    CDS = data['Mutation CDS']
    all_number = len(data.index)
    num = Num - 1
    unit = all_number // num
    max_index = index * unit
    start_index = max_index - unit
    if(index <= num):
        temp_CDS = CDS[start_index:max_index]
        temp_data = data[start_index:max_index]
    else:
        temp_CDS = CDS[start_index:]
        temp_data = data[start_index:]
    del CDS
    del data

    return temp_data


main()
