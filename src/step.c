#include <stdlib.h>
#include "step.h"
#include "hydro/hydro.h"

void substep_forward(double *prim1, double *prim2, double r1, double r2)
{
    int nc = numc();
    int i;

    double *dprim = (double *) malloc(nc * sizeof(double));

    flow_grad(prim1, r1, dprim);
    for(i=0; i<nc; i++)
        prim2[i] = prim1[i] + (r2-r1)*dprim[i];

    free(dprim);
}


