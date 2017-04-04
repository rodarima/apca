import numpy as np
from sys import stdin

data = np.genfromtxt(stdin.buffer, delimiter="\t")
#data = np.genfromtxt('firsthh.csv', delimiter='\t')

cols = data.shape[1]

if cols != 3:
	print("Data has {} columns, but I expected 3".format(cols))
	exit(1)

# First classify by the size of the matrix

size_set = np.unique(np.sort(data[:,1]))
print(size_set)
exit()

for size in size_set:
	var_data = data[data[:,1] == size]
#	print(var, var_data.shape)

	bit_col = np.unique(np.sort(var_data[:,0]))
#	print(bit_col)
	for bits in bit_col:
		chunk = var_data[var_data[:,0] == bits]
		y = chunk[:,2]
		mean = np.mean(y)
		print("{}\t{}\t{}".format(int(var), int(bits), mean))

exit()


# Sort by the bit-width

data_sort = data[data[:,0].argsort()]

x = np.unique(data_sort[:,0])
#x = np.sort(x)
#print(x)

#x = np.arange(2, 101, 5)
#print(x)

#runs = int(len(data[:,0])/len(np.unique(data_sort[:,0])))
#print(runs)

#mat = np.zeros([runs, len(x)])

for i in range(len(x)):
	bits = x[i]
	chunk = data_sort[data_sort[:,0] == bits]
	y = chunk[:,1]
	var = np.var(y)
	mean = np.mean(y)
	print("{}\t{}\t{}".format(bits, mean, var))
#	mat[:,i] = chunk[:,1]
#	median = np.median(y)
#	q25, q75 = np.percentile(y, [25, 75])
#	iqr = q75 - q25
#	lw = np.min(y[y > q25 - 1.5*iqr])
#	uw = np.max(y[y < q75 + 1.5*iqr])
#	print(median,q25,q75,lw,uw)
#	exit()

#np.savetxt("boxplot.csv", mat, delimiter=",")
