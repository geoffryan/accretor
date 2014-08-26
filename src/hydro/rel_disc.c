#include <math.h>
#include "hydro.h"
#include "../eos/eos.h"

int numq_rel_disc()
{
    return 6;
}

int numc_rel_disc()
{
    return 6;
}

void initial_rel_disc(double *prim, double *R1, double *R2)
{
    *R1 = 100.0;
    *R2 = 10.0;

    prim[RHO] = 1.0;
    prim[TTT] = 1.0;
    prim[URR] = -0.000001;
    prim[UPP] = 0.001;
    prim[DUR] = 0.0;
    prim[DUP] = 0.0;
}

void flow_grad_rel_disc(double *prim, double r, double *dprim)
{
    double drho, dttt, dur, dup, ddur, ddup;
    double u0, ur, up, u0d, urd, upd, du0, s00, s0r, s0p, srr, srp, spp, s2;
    double psi, dpsi, ddpsi, P, nu, eps, deps, h, qdot;

    psi = 1.0 - 2*M/r;
    dpsi = 2*M/(r*r);
    ddpsi = -4*M/(r*r*r);

    ur = prim[URR];
    up = prim[UPP];
    u0 =  sqrt((1.0 + ur*ur/psi + r*r*up*up)/psi);

    u0d = -psi*u0;
    urd = ur/psi;
    upd = r*r*up;

    dur = prim[DUR];
    dup = prim[DUP];
    du0 = 0.5*u0*((2*ur*dur/psi-ur*ur*dpsi/(psi*psi)+2*r*up*up+2*r*r*up*dup)*psi - (1.0+ur*ur/psi+r*r*up*up)*dpsi) / (psi*psi);

    srr = (ur*ur+psi)*dur - 2*r*psi*ur*up*up - (1/r+dpsi/psi)*ur*ur*ur + (u0*u0*psi*dpsi-dpsi-psi/r)*ur;
    srp = (ur*ur+psi)*dup - r*psi*up*up*up + 0.5*psi*dpsi*u0*u0 + (1/r-0.5*dpsi/psi)*ur*ur*up;
    spp = 2*ur*up*dup - (1/(r*r)+up*up)*dur + 3*ur*up*up/r+ur/(r*r*r);

    s0r = -(urd*srr+upd*srp)/u0d;
    s0p = -(urd*srp+upd*spp)/u0d;
    s00 = -(urd*s0r+upd*s0p)/u0d;

    s2 = (-psi)*(-psi)*s00*s00 + (1.0/psi)*(1.0/psi)*srr*srr + (r*r)*(r*r)*spp*spp
        +2*(-psi)*(1.0/psi)*s0r*s0r + 2*(-psi)*(r*r)*s0p*s0p + 2*(1.0/psi)*(r*r)*srp*srp;

    P = pressure(prim, r);
    eps = spec_int_en(prim, r);
    h = 1.0 + eps + P/prim[RHO];
    nu = alpha * sqrt(r*r*r/M) * P / (u0 * prim[RHO] * h);
    qdot = cool(prim, r);

    drho = -dur - ur/r;
    deps = (0.5*prim[RHO]*nu*s2 + P*ur*drho/prim[RHO] - qdot) / (prim[RHO]*ur);


    dprim[RHO] = drho;
    dprim[TTT] = dttt;
    dprim[URR] = dur;
    dprim[UPP] = dup;
    dprim[DUR] = ddur;
    dprim[DUP] = ddup;
}
