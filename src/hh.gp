set logscale y 2
set format y '%g'
set terminal png
set output 'hh.png'

plot 	'hh.csv' using 1:2 with lines, \
	'hh.csv' using 1:(2**(1-$1)) with lines, \
	'hh.csv' using 1:(200*2**(1-$1)) with lines
