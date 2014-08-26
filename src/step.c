#include <stdlib.h>
#include <stdio.h>
#include "step.h"
#include "hydro/hydro.h"

void step_setup(int choice)
{
    if(choice == 0)
    {
        step = &forward_euler;
    }
    else if(choice == 1)
    {
        step = &forward_rk2;
    }
    else if(choice == 2)
    {
        step = &forward_rk4;
    }
}

void evolve(double *prim, double R1, double R2, int n)
{
    int i, nc, nq;
    double r = R1;
    double dr = (R2-R1)/n;

    nc = numc();
    nq = numq();

    while((R1<R2 && r<R2) || (R2<R1 && r>R2))
    {
        printf("r: %.12g, dr: %.12g\n", r, dr);

        step(prim, r, &dr);
        r += dr;

        FILE *f = fopen("out.txt", "a");
        fprintf(f, "%.12g", r);
        for(i=0; i<nc; i++)
            fprintf(f, " %.12g", prim[i]);
        fprintf(f, "\n");
        fclose(f);
    }
}

void substep_forward(double *prim1, double *prim2, double r, double dr)
{
    int nc = numc();
    int nq = numq();
    int i;

    double *dprim = (double *) malloc(nc * sizeof(double));

    flow_grad(prim1, r, dprim);
    for(i=0; i<nc; i++)
        prim2[i] = prim1[i] + dr*dprim[i];
    for(i=nc; i<nq; i++)
        prim2[i] = prim1[i];

    free(dprim);
}

void forward_euler(double *prim, double r, double *dr)
{
    substep_forward(prim, prim, r, *dr);
}

void forward_rk2(double *prim, double r, double *dr)
{
    int i;
    int nq = numq();
    int nc = numc();
    double *prim1 = (double *) malloc(nq * sizeof(double));
    double *prim2 = (double *) malloc(nq * sizeof(double));

    substep_forward(prim, prim1, r, 0.5*(*dr));
    substep_forward(prim1, prim2, r+0.5*(*dr), *dr);

    for(i=0; i<nc; i++)
        prim[i] = prim2[i] - prim1[i] + prim[i];

    free(prim1);
    free(prim2);
}

void forward_rk4(double *prim, double r, double *dr)
{
    int i;
    int nq = numq();
    int nc = numc();
    double *prim1 = (double *) malloc(nq * sizeof(double));
    double *prim2 = (double *) malloc(nq * sizeof(double));
    double *prim3 = (double *) malloc(nq * sizeof(double));
    double *prim4 = (double *) malloc(nq * sizeof(double));

    substep_forward(prim, prim1, r, 0.5*(*dr));

    substep_forward(prim1, prim2, r+0.5*(*dr), 0.5*(*dr));
    for(i=0; i<nc; i++)
        prim2[i] += prim[i] - prim1[i];
    
    substep_forward(prim2, prim3, r+0.5*(*dr), *dr);
    for(i=0; i<nc; i++)
        prim3[i] += prim[i] - prim2[i];
    
    substep_forward(prim3, prim4, r+(*dr), 0.5*(*dr));
    for(i=0; i<nc; i++)
        prim4[i] += prim[i] - prim3[i];

    for(i=0; i<nc; i++)
        prim[i] = (prim1[i]+2*prim2[i]+prim3[i]+prim4[i]-2*prim[i])/3.0;

    free(prim1);
    free(prim2);
    free(prim3);
    free(prim4);
}
