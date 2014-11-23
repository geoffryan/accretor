#include <math.h>
#include "eos.h"
#include "../hydro/hydro.h"

double pressure_ideal(double *prim, double r)
{
    return prim[RHO]*prim[TTT];
}

double spec_int_en_ideal(double *prim, double r)
{
    return prim[TTT]/(GAMMA-1.0);
}

double depsdrho_ideal(double *prim, double r)
{
    return 0.0;
}

double depsdT_ideal(double *prim, double r)
{
    return 1.0/(GAMMA-1.0);
}

double dPdrho_ideal(double *prim, double r)
{
    return prim[TTT];
}

double dPdT_ideal(double *prim, double r)
{
    return prim[RHO];
}
