#include <math.h>
#include <stdio.h>
#include "hydro.h"
#include "../eos/eos.h"

// Newtonian Disc ala Shakura-Sunyaev (ie. only r-phi viscous term)

int numq_newt_disc_SS()
{
    return 6;
}

int numc_newt_disc_SS()
{
    return 4;
}

void initial_newt_disc_SS(double *prim, double *R1, double *R2)
{
    double r1 = 3.0;
    double r2 = 2.0;
    double Mdot = 1.0e-6;
    double rho = 1.0;
    double risco = 1.0e-6;

    *R1 = r1;
    *R2 = r2;

    prim[RHO] = rho;
    prim[TTT] = sqrt(M/r1*r1*r1) * Mdot/(3*M_PI*alpha*rho);
    prim[URR] = -Mdot/(2*M_PI*r1*rho);
    prim[UPP] = sqrt(M/(r1*r1*r1));
    prim[ACC] = Mdot;
    //prim[LLL] = risco*risco*Mdot/(2*M_PI)*sqrt(M/(risco*risco*risco));
    prim[LLL] = 0.0;
}

void flow_grad_newt_disc_SS(double *prim, double r, double *dprim)
{
    double rho, T, v, om, P, Mdot, L;
    double drho, dv, dom, dT;
    double eps, nu, qdot;
    double dPPPdRHO, dPPPdTTT, dEPSdRHO, dEPSdTTT;
    double D;

    v = prim[URR];
    om = prim[UPP];
    T = prim[TTT];
    Mdot = prim[ACC];
    L = prim[LLL];

    rho = -Mdot / (2*M_PI*r*v);
    prim[RHO] = rho;

    P = pressure(prim, r);
    eps = spec_int_en(prim, r);
    nu = alpha * sqrt(r*r*r/M) * P / prim[RHO];
    qdot = cool(prim, r);
    dPPPdRHO = dPdrho(prim, r);
    dPPPdTTT = dPdT(prim, r);
    dEPSdRHO = depsdrho(prim, r);
    dEPSdTTT = depsdT(prim, r);

    D = v*v - dPPPdRHO - (P/(rho*rho) - dEPSdRHO)*dPPPdTTT/dEPSdTTT;

    // RHO derivative
    drho = 0.0;

    //om derivative
    dom = (r*r*r*rho*v*om - L) / (r*r*r*rho*nu);

    // v derivative
    dv = dPPPdRHO/r + r*om*om - M/(r*r) - (dEPSdRHO/r - P/(r*rho*rho) 
            + 0.5*r*nu*dom*dom/(rho*v) - qdot/(rho*rho*v))*dPPPdTTT/dEPSdTTT;
    dv *= v/D;

    // T derivative
    dT = v*v*(rho*dEPSdRHO-P/rho)/r + (rho*dEPSdRHO-P/rho)*(r*om*om-M/(r*r))
        + (v*v-dPPPdRHO)*(0.5*r*r*nu*dom*dom-qdot/rho)/v;
    dT /= dEPSdTTT * D;

    printf("        %.10lg %.10lg %.10lg %.10lg\n", drho, dT, dv, dom);

    dprim[RHO] = drho;
    dprim[TTT] = dT;
    dprim[URR] = dv;
    dprim[UPP] = dom;
}
