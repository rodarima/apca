import pandas as p
import numpy as np
from sys import stdin,argv

degree = 1
if len(argv) == 2:
	degree = int(argv[1])

data = np.genfromtxt(stdin.buffer, delimiter=" ")
cols = data.shape[1]

if cols != 2:
	print("Data has {} columns, but I expected 2".format(cols))
	exit(1)

x = data[:,0]
y = data[:,1]

regression = np.polyfit(x, y, degree)

print(regression)
