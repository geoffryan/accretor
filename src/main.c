#include <stdlib.h>
#include <stdio.h>
#include "step.h"
#include "hydro/hydro.h"

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

