#! /bin/bash

F=${1%.csv}

awk '{print $1 "\t" $3 "\t" $6}' < ${F}.csv > ${F}a.csv
awk '{print $2 "\t" ($1+log($3)/log(2))}' < ${F}a.csv > ${F}b.csv
python py/stat.py < ${F}b.csv > ${F}c.csv
awk '{print $1 "\t" ($2*$2*$2)}' ${F}c.csv > ${F}d.csv
