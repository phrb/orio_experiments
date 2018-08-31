#! /bin/bash

for i in `seq 1 $1`;
do
	../../../orcc -v -e $2
done

mkdir xeon_e5_2630_v3_$(uname -n | cut -d. -f1)

./db2csv.py
mv results.* search_space* xeon_e5_2630_v3_$(uname -n | cut -d. -f1)

./clean.sh

git add --all
git commit -m "stencil3d Xeon E5 2630 v3 $(uname -n | cut -d. -f1) $(date)"
git push
