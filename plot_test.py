import sys
import matplotlib.pyplot as plt
import numpy as np

fname = "out.txt"

if len(sys.argv) == 2:
    fname = sys.argv[1]

x, y = np.loadtxt(fname, unpack=True, usecols=[0,1])

X = np.linspace(x.min(), x.max(), 500)
Y = np.cos(X)

plt.figure(figsize=(12,4))
plt.plot(X, Y, 'b')
plt.plot(x, y, 'k+')

plt.show()
