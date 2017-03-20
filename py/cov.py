import numpy as np
from scipy.stats.stats import pearsonr
from sys import stdin

data = np.genfromtxt(stdin.buffer, delimiter="\t")

cols = data.shape[1]

if cols != 2:
	print("Data has {} columns, but I expected 2".format(cols))
	exit(1)

x = data[:,0]
y = data[:,1]
print("Covariance matrix")
print(np.cov(x, y))
print("Correlation matrix")
print(np.corrcoef(x, y))
pcoef, pvalue = pearsonr(x, y)
print("Pearson coefficient: {}".format(pcoef))
print("p-value: {}".format(pvalue))

