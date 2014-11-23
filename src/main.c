#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "par.h"
#include "step.h"
#include "hydro/hydro.h"

int main(int argc, char *argv[])
{
    int i, nq, nc;
    double R1, R2;
    double *prim;

    struct parList theParList = PAR_DEFAULT;

    if(argc == 2)
        read_pars(&theParList, argv[1]);
    else
        read_pars(&theParList, "in.par");

    printf("%d\n", theParList.hydro);
    printf("%d\n", theParList.eos);
    printf("%d\n", theParList.cool);
    printf("%d\n", theParList.step);
    printf("%d\n", theParList.N);
    printf("%d\n", theParList.stop);

    printf("Setting up...\n");
    hydro_setup(&theParList);
    step_setup(&theParList);
    eos_setup(&theParList);
    cool_setup(&theParList);
    
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
    evolve(prim, R1, R2, &theParList);
    
    printf("Cleaning...\n");
    free(prim);
    return 0;
}

