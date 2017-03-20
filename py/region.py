import numpy as np
from scipy.stats import norm
from sys import stdin


data = np.genfromtxt(stdin.buffer, delimiter="\t")
#data = np.genfromtxt('firsthh.csv', delimiter='\t')

if len(data.shape) != 1:
	print("Data has more than 1 column")
	exit(1)

n = len(data)
mean = np.mean(data)
var = np.var(data)
ALPHAS = [0.10, 0.05, 0.01, 0.005, 0.001]
print("Mean:     {:.5f}".format(mean))
print("Variance: {:.5f}".format(var))
for alpha in ALPHAS:
	semi = norm.ppf(1 - alpha/2) * np.sqrt(var/n)
	x_min = mean - semi
	x_max = mean + semi

	#print("{:.3f} ({:.5f}, {:.5f})".format(alpha, x_min, x_max))
	print("{:.3f}\t{:.5f} ± {:.5f}".format(alpha, mean, semi))

