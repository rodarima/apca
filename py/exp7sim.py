import numpy as np
import datetime
from sys import argv

n = 100
b = 10
nodes = 24

#if len(argv) >= 2:
#	n = int(argv[1])
#
#if len(argv) >= 3:
#	b = int(argv[2])

# Coefficients for Gold
G = np.array([  -21.89924143, 3.06483927])

# Coefficients for b bits
B = np.array([  -21.44599281, 3.04707157])

#t = n**G[1] * 2**G[0] + b * n**G[1] * 2**G[0]
#print(t)
#print(str(datetime.timedelta(seconds=t)))

def time_sim(n, b):
	return n**G[1] * 2**G[0] + b * n**G[1] * 2**G[0]

N = range(600, 5000+1, 200)
r = 5

cpu_time = 0

for n in N:
	cpu_time += r * time_sim(n, b)

wall_time = cpu_time / nodes

print([i for i in N])
print(wall_time)
print(str(datetime.timedelta(seconds=wall_time)))
