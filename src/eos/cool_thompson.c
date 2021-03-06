#include "eos.h"
#include "../hydro/hydro.h"

double cool_thompson(double *prim, double r)
{
    return q0 * prim[TTT]*prim[TTT]*prim[TTT]*prim[TTT] / prim[RHO];
}
