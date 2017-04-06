import numpy as np
from sys import stdin

data = np.genfromtxt(stdin.buffer, delimiter=" ")
#data = np.genfromtxt('firsthh.csv', delimiter='\t')

cols = data.shape[1]

# n bits err time time_gold

if cols != 5:
	print("Data has {} columns, but I expected 5".format(cols))
	exit(1)

# First classify by the size n of the matrix

n_col = np.unique(np.sort(data[:,0]))
#print(var_col)

for n in n_col:
	n_data = data[data[:,0] == n]
#	print(var, var_data.shape)

	bit_col = np.unique(np.sort(n_data[:,1]))
#	print(bit_col)
	for bits in bit_col:
		chunk = n_data[n_data[:,1] == bits]
		err = chunk[:,2]
		time = chunk[:,3]
		time_gold = chunk[:,4]
		print("{} {} {} {}".format(
			int(n), int(bits), np.mean(err), np.mean(time)))

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
