import numpy as np
from sys import stdin

data = np.genfromtxt(stdin.buffer, delimiter=" ", skip_header=1)
#data = np.genfromtxt('firsthh.csv', delimiter='\t')

cols = data.shape[1]

if cols != 3:
	print("Data has {} columns, but I expected 3".format(cols))
	exit(1)

# First classify by the variable affected

var_col = np.unique(np.sort(data[:,0]))
#print(var_col)

for var in var_col:
	var_data = data[data[:,0] == var]
#	print(var, var_data.shape)

	y = var_data[:,2]
	mean = np.mean(y)

	print("{}\t{}".format(int(var), mean))

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
