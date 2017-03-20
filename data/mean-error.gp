set format y '%g'
set terminal png #size 2000,1000
set output '../doc/img/mean-error.png'

set grid xtics
set grid ytics


plot 	'exp1a.csv' using 1:2 with lines
