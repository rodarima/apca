import numpy as np
from sys import stdin, argv

data = np.genfromtxt(stdin.buffer, delimiter=" ")


rows, cols = data.shape

col = int(argv[1]) - 1
samples = 5
data_sort = data[data[:,col].argsort(), :]

dataset = data_sort

num = 0
last = None

for i in range(rows):
	row = dataset[i,:]
	if row[col] != last:
		num = 0
		last = row[col]
	elif num < samples:
		num += 1
		print(" ".join(map(str, row)))

