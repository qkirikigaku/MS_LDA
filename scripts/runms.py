import sys
import os
from multiprocessing import Pool
import subprocess

def main():
    arguments = []
    for i in range(1,51):
        for j in range(1,30):
            arguments.append((str(j), data_type, threshold, cancer_type, 
                              str(i)))
    pool = Pool()
    _ = pool.starmap(execute, arguments)

def execute(num_topic, data_type, threshold, cancer_type, experiment):
    result_path = 'result/data' + data_type + '_o' + threshold + '_' +\
                  cancer_type + '_' + experiment
    if(os.path.exists(result_path) == False): os.mkdir(result_path)
    cmd = 'bin/MS ' + num_topic + ' ' + data_type + ' ' + threshold +\
          ' ' + cancer_type + ' ' + experiment
    subprocess.call(cmd.split())

if __name__ == '__main__':
    args = sys.argv
    data_type = args[1]
    threshold = args[2]
    cancer_type = args[3]
    main()
