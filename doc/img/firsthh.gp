set logscale y 2
set format y '%g'
set terminal png size 2000,1000
set output 'firsthh.png'

set grid xtics
set grid ytics

set xrange [1:101]

set yrange [1e-34:1e2]

set xtics 2,4,100
set ytics 1e-34,100,100

plot 	'firsthh.csv' using 1:2

#plot 	'firsthh.csv' using 1:2, \
#	'firsthh.csv' using ($1+0.5):3, \
#	'firsthh.csv' using 1:(200*2**(1-$1)) with lines
