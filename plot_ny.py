import sys
import math
import matplotlib.pyplot as plt
import numpy as np

fname = "out.txt"

if len(sys.argv) == 2:
    fname = sys.argv[1]

r, rho, T, vr, vp = np.loadtxt(fname, unpack=True, usecols=[0,1,2,3,4])
f = open(fname, 'r')
line = f.readline()
md = float(line.split()[5])
f.close()

i = r.argmax()
R = np.linspace(r.min(), r.max(), 500)
RHO = rho[i] * np.power(R/r[i], -0.5)
TTT = T[i] * np.power(R/r[i], -1.5)
URR = vr[i] * np.power(R/r[i], -0.5)
UPP = vp[i] * np.power(R/r[i], -1.5)

xscale = 'log'
yscale = 'log'

plt.figure(figsize=(12,9))

plt.subplot(3,2,1)
plt.plot(r[:-1], rho[1:], 'k+')
plt.plot(R, RHO, 'r')
plt.xlabel(r"$r$")
plt.ylabel(r"$\rho$")
plt.gca().set_xscale(xscale)
plt.gca().set_yscale(yscale)

plt.subplot(3,2,2)
plt.plot(r, T, 'k+')
plt.plot(R, TTT, 'r')
plt.xlabel(r"$r$")
plt.ylabel(r"$T$")
plt.gca().set_xscale(xscale)
plt.gca().set_yscale(yscale)

plt.subplot(3,2,3)
if xscale == 'log':
    plt.plot(r, -vr, 'k+')
    plt.plot(R, -URR, 'r')
else:
    plt.plot(r, vr, 'k+')
    plt.plot(R, URR, 'r')
plt.xlabel(r"$r$")
plt.ylabel(r"$v_r$")
plt.gca().set_xscale(xscale)
plt.gca().set_yscale(yscale)

plt.subplot(3,2,4)
plt.plot(r, vp, 'k+')
plt.plot(R, UPP, 'r')
plt.xlabel(r"$r$")
plt.ylabel(r"$v_\phi$")
plt.gca().set_xscale(xscale)
plt.gca().set_yscale(yscale)

plt.subplot(3,2,5)
plt.hlines(md,r[0], r[-1], color='k', lw=2)
#plt.plot(r[:-1], -2*math.pi*r[:-1]*rho[1:]*vr[:-1], 'b^', mew=0)
plt.plot(r, -2*math.pi*r*rho*vr, 'b^', mew=0)
plt.plot(R, -2*math.pi*R*RHO*URR, 'r')
plt.xlabel(r"$r$")
plt.ylabel(r"$\dot{\mathcal{M}}$")
#plt.gca().set_xscale(xscale)
#plt.gca().set_yscale(yscale)
plt.tight_layout()

plt.show()
