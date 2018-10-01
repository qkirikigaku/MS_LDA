#!/bin/bash
division=$1
echo "make d3..."
for i in $(seq 1 ${division}); do
    python Preprocessing/make_d3.py ${i} ${division}
done
python Preprocessing/integrate_d3.py 400 ${division} 
