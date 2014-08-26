#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "step.h"
#include "hydro/hydro.h"

int main(int argc, char *argv[])
{
    int i, nq, nc, N;
    double R1, R2;
    double *prim;

    N = 128;

    printf("Setting up...\n");
    hydro_setup(0);
    step_setup(2);
    
    printf("Initializing...\n");
    nq = numq();
    nc = numc();
    prim = (double *) malloc(nq * sizeof(double));
    initial(prim, &R1, &R2);

    printf("Printing...\n");
    FILE *f = fopen("out.txt", "w");
    fprintf(f, "%.12g", R1);
    for(i=0; i<nq; i++)
        fprintf(f, " %.12g", prim[i]);
    fprintf(f, "\n");
    fclose(f);

    printf("Evolving...\n");
    evolve(prim, R1, R2, N);
    
    printf("Cleaning...\n");
    free(prim);
    return 0;
}

