#include <math.h>
#include <stdio.h>
#include "hydro.h"
#include "../eos/eos.h"

// Disc in the Schwarzschild metric.  Approximating the Popham-Gammie 98
// setup.

int numq_rel_disc_PG()
{
    return 6;
}

int numc_rel_disc_PG()
{
    return 4;
}

void initial_rel_disc_PG(double *prim, double *R1, double *R2)
{
    *R1 = 1.0;
    *R2 = 3.0;

    prim[RHO] = 1.0e4;
    prim[TTT] = 1.0e-8;
    prim[URR] = -1.0e-7;
    prim[UPP] = 0.001;
    prim[DUR] = 1.0e-7;
    prim[DUP] = -0.0015;
}

void flow_grad_rel_disc_PG(double *prim, double r, double *dprim)
{

    double drho, dttt, dur, dup, ddur, ddup;
    double u0, ur, up, u0d, urd, upd, du0, s00, s0r, s0p, srr, srp, spp, s2, dsrp, dsrr;
    double psi, dpsi, ddpsi, grr, dgrr;
    double P, dP, nu, dnu, eps, deps, h, dh, qdot;

    psi = 1.0 - 2*M/r;
    dpsi = 2*M/(r*r);
    ddpsi = -4*M/(r*r*r);
    grr = 1.0/psi;
    dgrr = -dpsi/(psi*psi);

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

    // RHO derivative
    drho = -prim[RHO]*(dur + ur/r);
   
   // printf("s00: %.12g srr: %.12g spp: %.12g\n", s00, srr, spp);
   // printf("s0r: %.12g s0p: %.12g srp: %.12g\n", s0r, s0p, srp);

    // TTT derivative
    deps = (0.5*prim[RHO]*nu*s2 + P*ur*drho/(prim[RHO]*prim[RHO]) - qdot) / (prim[RHO]*ur);
    dttt = (deps - depsdrho(prim,r)*drho) / depsdT(prim,r); 

    // Derivatives of P,nu,h
    dP = dPdrho(prim,r)*drho + dPdT(prim,r)*dttt;
    dh = deps + dP/prim[RHO] - P*drho/(prim[RHO]*prim[RHO]);
    dnu = nu*(-du0/u0 + 1.5/r + dP/P - (prim[RHO]*dh+drho*h)/(prim[RHO]*h));

    // Derivative of r-phi stress tensor component
    dsrp = (prim[RHO]*ur*(h*(2*r*up+r*r*dup)+upd*dh) + upd*qdot - 3*r*srp*prim[RHO]*nu - r*r*srp*(drho*nu+prim[RHO]*dnu)) / (r*r*prim[RHO]*nu);

    // DUP derivative (from Mathematica calculation)
    ddup = dsrp;
    ddup -= psi*dpsi*u0*du0*up - (psi+r*dpsi)*up*up*up + (2/r-dpsi/psi)*ur*up*dur
         + (-3*r*up*up*psi + 2*ur*dur + dpsi + 0.5*u0*u0*psi*dpsi + ur*ur*(1/r-0.5*dpsi/psi))*dup
          + up*ur*ur*(-1.0/(r*r) + 0.5*dpsi*dpsi/(psi*psi) - 0.5*ddpsi/psi)
          + 0.5*up*u0*u0*(dpsi*dpsi + psi*ddpsi);
    ddup /= ur*ur+psi;


    // Derivative of r-r stress tensor component
    dsrr = prim[RHO]*ur*(h*ur*dgrr+h*dur*grr+dh*ur*grr);
    dsrr += urd*qdot;
    dsrr -= prim[RHO]*nu*grr*srr/r + grr*srr*(prim[RHO]*dnu+drho*nu) + prim[RHO]*nu*srr*dgrr;
    dsrr += dP - 0.5*prim[RHO]*((h*u0*u0-nu*s00)*(-dpsi)
            + (h*ur*ur-nu*srr)*dgrr + (h*up*up-nu*spp)*(2*r));
    dsrr /= prim[RHO]*nu*grr;

    // DUR Derivative (from Mathemtica calculation)
    ddur = dsrr;
    ddur -= 2*ur*dur*dur - 4*r*psi*ur*up*dup + 2*psi*dpsi*u0*ur*du0
            - 2*ur*up*up*(psi+r*dpsi) + ur*ur*ur*(1.0/(r*r)+dpsi*dpsi/(psi*psi)-ddpsi/psi)
            + (-psi/r - 2*r*psi*up*up + psi*dpsi*u0*u0 - 3*ur*ur/r - 3*ur*ur*dpsi/psi)*dur
            + ur*(psi/(r*r) - dpsi/r - ddpsi + u0*u0*(dpsi*dpsi+psi*ddpsi));
    ddur /= ur*ur + psi;

    dprim[RHO] = drho;
    dprim[TTT] = dttt;
    dprim[URR] = dur;
    dprim[UPP] = dup;
    dprim[DUR] = ddur;
    dprim[DUP] = ddup;
}
