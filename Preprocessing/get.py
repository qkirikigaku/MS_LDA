import subprocess


def main():
    cancer_types = load_cancer_type()
    for cancer_type in cancer_types:
        for dic in range(1, 5):
            cmd = 'python Preprocessing/get_M' + str(dic) + '.py 400 ' + \
                cancer_type
            subprocess.call(cmd.split())


def load_cancer_type():
    file_name = 'data/cancer_types.txt'
    lines = open(file_name, 'r').readlines()
    cancer_types = list()
    for line in lines:
        cancer_types.append(line)
    return cancer_types


if __name__ == '__main__':
    main()
