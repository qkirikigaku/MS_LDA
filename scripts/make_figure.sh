#!/bin/bash
data_type=$1
threshold=$2
cancer_type=$3
str1=matsutani@133.9.8.88:project/MS_LDA/result/data${data_type}_o${threshold}_${cancer_type}*
scp -r ${str1} result
python Drawing/find_best_K.py ${data_type} ${threshold} ${cancer_type}
FILENAME=ref/${data_type}_${threshold}_${cancer_type}.txt
cnt=0
array=()
while read line;
do
    cnt=$(expr $cnt + 1)
    if test $cnt -eq 1; then
	num_data=$line
    fi
    if test $cnt -eq 2; then
	number_of_topic=$line
    fi
    if test $cnt -ge 3; then
	array+=($line)
    fi
done<$FILENAME
echo ${number_of_topic}
for i in $(seq 2 30); do
    d_num=${array[i-2]}
    str1=result/data${data_type}_o${threshold}_${cancer_type}_${d_num}
    str2=../data${data_type}_o${threshold}_${cancer_type}_${num_data}
    cd ${str1}
    if test ${i} -le 9; then
	str3=result_k0${i}.txt
    else
	str3=result_k${i}.txt
    fi
    cp ${str3} ${str2}
    cd ../..
done
for i in $(seq 1 50); do
    str1=result/data${data_type}_o${threshold}_${cancer_type}_${i}
    if test ${i} -eq ${num_data}; then
	cd ${str1}
	mkdir figure
	cd ../..
	python Drawing/comparison_K.py ${data_type} ${threshold} ${cancer_type} ${i}
	str2=result/data${data_type}_o${threshold}_${cancer_type}
	mv ${str1} ${str2}
    else
	rm -r ${str1}
    fi
done
str2=result/data${data_type}_o${threshold}_${cancer_type}
cd ${str2}
cd figure
mkdir ${number_of_topic}
cd ../../..
sh scripts/make_draw.sh ${data_type} ${number_of_topic} \
    ${threshold} ${cancer_type}
