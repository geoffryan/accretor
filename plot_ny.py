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
R = np.linspace(r.min(), r.max(), 2000)
RHO = rho[i] * np.power(R/r[i], -0.5)
TTT = T[i] * np.power(R/r[i], -1.0)
URR = vr[i] * np.power(R/r[i], -0.5)
UPP = vp[i] * np.power(R/r[i], -1.5)

xscale = 'linear'
yscale = 'linear'

alpha = 0.01
M = 1.0
GAMMA = 4.0/3.0
f = 1.0

nu = alpha * np.sqrt(r*r*r/M) * T
NU = alpha * np.sqrt(R*R*R/M) * TTT

frho = r[1:-1]*rho[1:-1]*vr[1:-1]
fsr = r[1:-1]*rho[1:-1]*vr[1:-1]*vr[1:-1] + r[1:-1]*rho[1:-1]*T[1:-1]
fsp = r[1:-1]*r[1:-1]*r[1:-1]*(rho[1:-1]*vr[1:-1]*vp[1:-1]
        - rho[1:-1]*nu[1:-1]*(vp[2:]-vp[:-2])/(r[2:]-r[:-2]))
fe = r[1:-1]*rho[1:-1]*T[1:-1]*vr[1:-1] / (GAMMA-1.0)

FRHO = R[1:-1]*RHO[1:-1]*URR[1:-1]
FSRR = R[1:-1]*RHO[1:-1]*URR[1:-1]*URR[1:-1] + R[1:-1]*RHO[1:-1]*TTT[1:-1]
FSPP = R[1:-1]*R[1:-1]*(R[1:-1]*RHO[1:-1]*URR[1:-1]*UPP[1:-1]
        + 1.5*RHO[1:-1]*NU[1:-1]*UPP[1:-1])
FTAU = R[1:-1]*RHO[1:-1]*TTT[1:-1]*URR[1:-1] / (GAMMA-1.0)

srho = np.zeros(len(frho)-2)
ssr = r[2:-2]*rho[2:-2]*vp[2:-2]*vp[2:-2] + rho[2:-2]*T[2:-2]/r[2:-2] - M*rho[2:-2] / (r[2:-2]*r[2:-2])
ssp = np.zeros(len(frho)-2)
se = -rho[2:-2]*T[2:-2]*((r[3:-1]*vr[3:-1]-r[1:-3]*vr[1:-3])/(r[3:-1]-r[1:-3]))/r[2:-2] \
        + 0.5*f*rho[2:-2]*nu[2:-2]*r[2:-2]*r[2:-2] \
        * 2*((vp[3:-1]-vp[1:-3])**2)/((r[3:-1]-r[1:-3])**2)

SRHO = np.zeros(len(FRHO)-2)
SSRR = R[2:-2]*RHO[2:-2]*UPP[2:-2]*UPP[2:-2] + RHO[2:-2]*TTT[2:-2]/R[2:-2] - M*RHO[2:-2] / (R[2:-2]*R[2:-2])
SSPP = np.zeros(len(FRHO)-2)
STAU = -RHO[2:-2]*TTT[2:-2]*((R[3:-1]*URR[3:-1]-R[1:-3]*URR[1:-3])/(R[3:-1]-R[1:-3]))/R[2:-2] \
        + 0.5*f*RHO[2:-2]*NU[2:-2]*R[2:-2]*R[2:-2] \
        * 2*((UPP[3:-1]-UPP[1:-3])**2)/((R[3:-1]-R[1:-3])**2)

fig1 = plt.figure(figsize=(12,9))

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

fig2 = plt.figure(figsize=(12,9))

plt.subplot(2,2,1)
plt.plot(r[2:-2], (frho[2:]-frho[:-2])/(r[3:-1]-r[1:-3]) / r[2:-2], 'k,')
plt.plot(r[2:-2], srho, 'b,')
plt.plot(R[2:-2], (FRHO[2:]-FRHO[:-2])/(R[3:-1]-R[1:-3]) / R[2:-2], 'g,')
plt.plot(R[2:-2], SRHO, 'r,')
plt.xlabel(r"$r$")
plt.ylabel(r"$F_\rho$")
plt.gca().set_xscale(xscale)
plt.gca().set_yscale(yscale)

plt.subplot(2,2,2)
plt.plot(r[2:-2], (fsr[2:]-fsr[:-2])/(r[3:-1]-r[1:-3]) / r[2:-2], 'k,')
plt.plot(r[2:-2], ssr, 'b,')
plt.plot(R[2:-2], (FSRR[2:]-FSRR[:-2])/(R[3:-1]-R[1:-3]) / R[2:-2], 'g,')
plt.plot(R[2:-2], SSRR, 'r,')
plt.xlabel(r"$r$")
plt.ylabel(r"$F_{Sr}$")
plt.gca().set_xscale(xscale)
plt.gca().set_yscale(yscale)

plt.subplot(2,2,3)
plt.plot(r[2:-2], (fsp[2:]-fsp[:-2])/(r[3:-1]-r[1:-3]) / r[2:-2], 'k,')
plt.plot(r[2:-2], ssp, 'b,')
plt.plot(R[2:-2], (FSPP[2:]-FSPP[:-2])/(R[3:-1]-R[1:-3]) / R[2:-2], 'g,')
plt.plot(R[2:-2], SSPP, 'r,')
plt.xlabel(r"$r$")
plt.ylabel(r"$F_{S\phi}$")
plt.gca().set_xscale(xscale)
plt.gca().set_yscale(yscale)

plt.subplot(2,2,4)
plt.plot(r[2:-2], (fe[2:]-fe[:-2])/(r[3:-1]-r[1:-3]) / r[2:-2], 'k,')
plt.plot(r[2:-2], se, 'b,')
plt.plot(R[2:-2], (FTAU[2:]-FTAU[:-2])/(R[3:-1]-R[1:-3]) / R[2:-2], 'g,')
plt.plot(R[2:-2], STAU, 'r,')
plt.xlabel(r"$r$")
plt.ylabel(r"$F_{\tau}$")
plt.gca().set_xscale(xscale)
plt.gca().set_yscale(yscale)

plt.tight_layout()

plt.show()
