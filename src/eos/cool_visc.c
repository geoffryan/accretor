#include <math.h>
#include "eos.h"
#include "../hydro/hydro.h"

double f = 0.5;

double cool_visc(double *prim, double r)
{
    return 2.25 * f * alpha * sqrt(r*r*r/M) * prim[RHO]*prim[TTT] * prim[UPP]*prim[UPP];
}
