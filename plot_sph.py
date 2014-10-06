import sys
import matplotlib.pyplot as plt
import numpy as np

fname = "out.txt"

if len(sys.argv) == 2:
    fname = sys.argv[1]

r, rho, T, vr = np.loadtxt(fname, unpack=True, usecols=[0,1,2,3])

GAM = 4.0/3.0
i = r.argmax()
R = np.linspace(r.min(), r.max(), 500)
RHO = rho[i] * np.power(R/r[i], -0.5)
TTT = T[i] * np.power(R/r[i], -1.0)
URR = vr[i] * np.power(R/r[i], -0.5)

xscale = 'log'
yscale = 'log'

plt.figure(figsize=(12,9))

plt.subplot(2,2,1)
plt.plot(r, rho, 'k')
#plt.plot(R, RHO, 'r')
plt.gca().set_xscale(xscale)
plt.gca().set_yscale(yscale)
plt.xlabel(r"$r$")
plt.ylabel(r"$\rho$")

plt.subplot(2,2,2)
plt.plot(r, T, 'k')
#plt.plot(R, TTT, 'r')
plt.gca().set_xscale(xscale)
plt.gca().set_yscale(yscale)
plt.xlabel(r"$r$")
plt.ylabel(r"$T$")

plt.subplot(2,2,3)
plt.plot(r, -vr, 'k')
plt.plot(r, np.sqrt(GAM*T), 'b')
#plt.plot(R, URR, 'r')
plt.gca().set_xscale(xscale)
plt.gca().set_yscale(yscale)
plt.xlabel(r"$r$")
plt.ylabel(r"$-v_r$")

plt.tight_layout()

plt.show()
