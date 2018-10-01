def main():
    for i in range(22):
        chr_file = 'raw_data/chroms/chr' + str(i+1) + '.fa'
        chr = open(chr_file, 'r')
        line1 = chr.readlines()
        flag = 0
        output_file = 'raw_data/chr' + str(i+1) + '.fa'
        output = open(output_file, 'w')
        for line in line1:
            if(flag == 0):
                flag = 1
            else:
                content = line.upper()
                output.write(content)
        chr.close()
        output.close()

    chr_file = 'raw_data/chroms/chrM.fa'
    chr= open(chr_file, 'r')
    line1 = chr.readlines()
    flag = 0
    output_file= 'raw_data/chrM.fa'
    output = open(output_file, 'w')
    for line in line1:
        if(flag == 0):
            flag = 1
        else:
            content = line.upper()
            output.write(content)
    chr.close()
    output.close()

    chr_file = 'raw_data/chroms/chrX.fa'
    chr= open(chr_file, 'r')
    line1 = chr.readlines()
    flag = 0
    output_file= 'raw_data/chrX.fa'
    output = open(output_file, 'w')
    for line in line1:
        if(flag == 0):
            flag = 1
        else:
            content = line.upper()
            output.write(content)
    chr.close()
    output.close()

    chr_file = 'raw_data/chroms/chrY.fa'
    chr= open(chr_file, 'r')
    line1 = chr.readlines()
    flag = 0
    output_file= 'raw_data/chrY.fa'
    output = open(output_file, 'w')
    for line in line1:
        if(flag == 0):
            flag = 1
        else:
            content = line.upper()
            output.write(content)
    chr.close()
    output.close()

if __name__ == '__main__':
    main()
