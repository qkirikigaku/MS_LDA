#!/bin/bash
data_type=$1
number_of_topic=$2
threshold=$3
cancer_type=$4

if [ ${data_type} -eq 1 ]; then
    python Drawing/matching.py ${data_type} ${number_of_topic} \
        ${threshold} ${cancer_type}
    python Drawing/data_represent.py ${data_type} ${number_of_topic} \
        ${threshold} ${cancer_type}
fi

if [ ${data_type} -eq 2 ]; then
    python Drawing/n2_to_n1.py ${threshold} ${number_of_topic} ${cancer_type}
    python Drawing/matching.py ${data_type} ${number_of_topic} \
        ${threshold} ${cancer_type}
    python Drawing/data_represent.py ${data_type} ${number_of_topic} \
        ${threshold} ${cancer_type}
    python Drawing/draw_data2.py ${data_type} ${number_of_topic} \
        ${threshold} ${cancer_type}
fi

if [ ${data_type} -eq 3 ]; then
    python Drawing/n3_to_n1.py ${threshold} ${number_of_topic} ${cancer_type}
    python Drawing/matching.py ${data_type} ${number_of_topic} \
        ${threshold} ${cancer_type}
    python Drawing/data_represent.py ${data_type} ${number_of_topic} \
        ${threshold} ${cancer_type}
    python Drawing/draw_data3_indel.py ${data_type} ${number_of_topic} \
        ${threshold} ${cancer_type}
fi
if [ ${data_type} -eq 4 ]; then
    python Drawing/n4_to_n1.py ${threshold} ${number_of_topic} ${cancer_type}
    python Drawing/matching.py ${data_type} ${number_of_topic} \
        ${threshold} ${cancer_type}
    python Drawing/data_represent.py ${data_type} ${number_of_topic} \
        ${threshold} ${cancer_type}
    python Drawing/draw_data4_indel.py ${data_type} ${number_of_topic} \
        ${threshold} ${cancer_type}
    python Drawing/draw_data4.py ${data_type} ${number_of_topic} \
        ${threshold} ${cancer_type}
fi

