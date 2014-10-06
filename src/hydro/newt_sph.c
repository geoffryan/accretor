#include <math.h>
#include <stdio.h>
#include "hydro.h"
#include "../eos/eos.h"

// Newtonian inviscid spherical flow

int numq_newt_sph()
{
    return 4;
}

int numc_newt_sph()
{
    return 3;
}

void initial_newt_sph(double *prim, double *R1, double *R2)
{
    double r1 = 10.0;
    double r2 = 0.01;
    double Mdot = 1.825*M_PI;
    double rho, T, vr;

    rho = 1.0;
    T = 1.0;
    vr = -Mdot/(4*M_PI*r1*r1*rho);

    *R1 = r1;
    *R2 = r2;

    prim[RHO] = rho;
    prim[TTT] = T;
    prim[URR] = vr;
    prim[SAC] = Mdot;
}

void exact_newt_sph(double *prim, double r)
{
    prim[RHO] = -prim[SAC] / (4*M_PI*r*r*prim[URR]);
}

void flow_grad_newt_sph(double *prim, double r, double *dprim)
{
    double rho, T, v, Mdot;
    double drho, dv, dT;
    double P, qdot;
    double dPPPdRHO, dPPPdTTT, dEPSdRHO, dEPSdTTT;
    double D;

    rho = prim[RHO];
    v = prim[URR];
    T = prim[TTT];
    Mdot = prim[SAC];

    P = pressure(prim, r);
    qdot = cool(prim, r);
    dPPPdRHO = dPdrho(prim, r);
    dPPPdTTT = dPdT(prim, r);
    dEPSdRHO = depsdrho(prim, r);
    dEPSdTTT = depsdT(prim, r);

    D = v*v - dPPPdRHO + (dEPSdRHO-P/(rho*rho))*dPPPdTTT/dEPSdTTT;

    // RHO derivative
    drho = 0.0;

    // v derivative
    dv = 2*(dPPPdRHO - (dEPSdRHO-P/(rho*rho))*dPPPdTTT/dEPSdTTT)*v/r
     - M*v/(r*r) + qdot*dPPPdTTT/(dEPSdTTT*rho*rho);
    dv /= D;

    // T derivative
    dT = (2*v*v/r-M/(r*r)) * (rho*dEPSdRHO-P/rho) - (v*v-dPPPdRHO)*qdot/(rho*v);
    dT /= dEPSdTTT * D;

 //    printf("        %.10lg %.10lg %.10lg\n", drho, dT, dv);

    dprim[RHO] = drho;
    dprim[TTT] = dT;
    dprim[URR] = dv;
}
