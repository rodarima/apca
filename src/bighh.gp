set logscale y 2
set format y '%g'
set terminal png size 2000,1000
set output 'bighh.png'

plot 	'bighh.csv' using 1:2, \
	'bighh.csv' using 1:(200*2**(1-$1)) with lines
