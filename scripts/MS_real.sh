#!/bin/bash
data_type=${1}
threshold=${2}
cancer_type=${3}
for i in $(seq 1 50); do
    python runms.py $i ${data_type} ${threshold} ${cancer_type}
done
