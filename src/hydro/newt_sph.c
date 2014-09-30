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
    double r1 = 3.0;
    double r2 = 1.0;
    double Mdot = 1.0e-6;
    double rho, T, vr;

    rho = 1.0;
    T = 1.0;
    vr = -Mdot/(2*M_PI*r1*rho);

    *R1 = r1;
    *R2 = r2;

    prim[RHO] = rho;
    prim[TTT] = T;
    prim[URR] = vr;
    prim[SAC] = Mdot;
}

void flow_grad_newt_sph(double *prim, double r, double *dprim)
{
    double rho, T, v, Mdot;
    double drho, dv, dT;
    double eps, qdot;
    double dPPPdRHO, dPPPdTTT, dEPSdRHO, dEPSdTTT;
    double D;

    v = prim[URR];
    T = prim[TTT];
    Mdot = prim[SAC];

    rho = -Mdot / (4*M_PI*r*r*v);
    prim[RHO] = rho;

    qdot = cool(prim, r);
    dPPPdRHO = dPdrho(prim, r);
    dPPPdTTT = dPdT(prim, r);
    dEPSdRHO = depsdrho(prim, r);
    dEPSdTTT = depsdT(prim, r);

    D = v*v - dPPPdRHO + dEPSdRHO*dPPPdTTT/dEPSdTTT;

    // RHO derivative
    drho = 0.0;

    // v derivative
    dv = 2*dPPPdRHO*v/r - M*v/(r*r) - (2*dEPSdRHO*v/r - qdot/(rho*rho))*dPPPdTTT/dEPSdTTT;
    dv /= D;

    // T derivative
    dT = 2*rho*v*v*dEPSdRHO/r - rho*dEPSdRHO*M/(r*r) -qdot*(v*v-dPPPdRHO)/(rho*v);
    dT /= dEPSdTTT * D;

 //    printf("        %.10lg %.10lg %.10lg\n", drho, dT, dv);

    dprim[RHO] = drho;
    dprim[TTT] = dT;
    dprim[URR] = dv;
}
