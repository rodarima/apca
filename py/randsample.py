import numpy as np
from sys import stdin, argv

data = np.genfromtxt(stdin.buffer, delimiter=" ")


rows, cols = data.shape

col = int(argv[1]) - 1
samples = 5
if len(argv) >= 3: samples = int(argv[2])

col_values = data[data[:,col].argsort(), col]
col_values = np.unique(col_values)

dataset = data

for val in col_values:
	possible_rows = dataset[dataset[:,col] == val, :]
	n = possible_rows.shape[0]
	indexes = np.arange(n)
	np.random.shuffle(indexes)
	for i in range(samples):
		row = possible_rows[indexes[i]]
		print(" ".join(map(str, row)))
