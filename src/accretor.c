#include <stdlib.h>
#include <stdio.h>
#include "accretor.h"
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

int main(int argc, char *argv[])
{
    int nq, nc;
    double *prim;

    printf("Setting up...\n");
    hydro_setup(0);
    
    printf("Initializing...\n");
    nq = numq();
    nc = numc();
    prim = (double *) malloc(nq * sizeof(double));
    initial(prim);

    printf("Printing...\n");
    int i;
    for(i=0; i<nq; i++)
        printf("%g\n", prim[i]);

    printf("Cleaning...\n");
    free(prim);
    return 0;
}

