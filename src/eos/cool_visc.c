#include <math.h>
#include "eos.h"
#include "../hydro/hydro.h"

double cool_visc(double *prim, double r)
{
	double rho, T, v, om, Mdot, L;
    double P, nu, dom;
    
    rho = prim[RHO];
    v = prim[URR];
    om = prim[UPP];
    T = prim[TTT];
    L = prim[LLL];
    P = pressure(prim, r);
    nu = alpha * sqrt(r*r*r/M) * P / prim[RHO];
	dom = (r*r*r*rho*v*om - L) / (r*r*r*rho*nu);

    return 0.5 * (1.0-visc_f) * alpha * sqrt(r*r*r/M) * rho*T * dom*dom;
}
