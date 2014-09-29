#include <math.h>
#include "eos.h"
#include "../hydro/hydro.h"

double GAM = 4.0/3.0;

double pressure_ideal(double *prim, double r)
{
    return prim[RHO]*prim[TTT];
}

double spec_int_en_ideal(double *prim, double r)
{
    return prim[TTT]/(GAM-1.0);
}

double depsdrho_ideal(double *prim, double r)
{
    return 0.0;
}

double depsdT_ideal(double *prim, double r)
{
    return 1.0/(GAM-1.0);
}

double dPdrho_ideal(double *prim, double r)
{
    return prim[TTT];
}

double dPdT_ideal(double *prim, double r)
{
    return prim[RHO];
}
