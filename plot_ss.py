import sys
import matplotlib.pyplot as plt
import numpy as np

fname = "out.txt"

if len(sys.argv) == 2:
    fname = sys.argv[1]

r, rho, T, vr, vp = np.loadtxt(fname, unpack=True, usecols=[0,1,2,3,4])

plt.figure(figsize=(12,9))

plt.subplot(2,2,1)
plt.plot(r, rho, 'k+')
plt.xlabel(r"$r$")
plt.ylabel(r"$\rho$")

plt.subplot(2,2,2)
plt.plot(r, T, 'k+')
plt.xlabel(r"$r$")
plt.ylabel(r"$T$")

plt.subplot(2,2,3)
plt.plot(r, vr, 'k+')
plt.xlabel(r"$r$")
plt.ylabel(r"$v_r$")

plt.subplot(2,2,4)
plt.plot(r, vp, 'k+')
plt.xlabel(r"$r$")
plt.ylabel(r"$v_\phi$")

plt.tight_layout()

plt.show()
