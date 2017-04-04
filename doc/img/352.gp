set logscale y 2
set format y '%g'
set terminal png
set output '352.png'

plot 	'352.csv' using 1:2 with lines, \
	'352.csv' using 1:(2**(1-$1)) with lines
