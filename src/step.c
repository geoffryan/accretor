#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "step.h"
#include "hydro/hydro.h"

/*
This file includes all of the integration algorithms for solving
a system of ODEs.  The current integrator is referred to by 
the "step" function pointer.  All integration algorithms have the
signature "void (*step)(double *prim, double r, double *dr)"
defined in step.h.  'prim' contains the current values of the 
dependent variables at location 'r', 'dr' refers to the desired step
size.  Adaptive integrators may alter the value of dr.  In calling
'step', the values in 'prim' will be updated to the values at r+dr.
*/

void step_setup(int choice)
{
    /*
    Intializes the 'step' function pointer to refer to the desired
    integration scheme.
    */

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
    else if(choice == 3)
    {
        step = &backward_euler;
    }
}

void evolve(double *prim, double R1, double R2, int n, int stop)
{
    /*
    Integrates prim from R1 to R2 in n (approximately) steps.
    */

    int i, iter, nc, nq;
    double r = R1;
    double dr = (R2-R1)/n;

    nc = numc();
    nq = numq();

    if(stop < 0)
        iter = -10;
    else
        iter = 0;

    while(((R1<R2 && r<R2) || (R2<R1 && r>R2)) && iter < stop)
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

        if(stop >= 0)
            iter++;
    }
}

void substep_forward(double *prim1, double *prim2, double r, double dr)
{
    /*
    Primitive explicit first-order step.  Uses prim1 as input, saves result
    in prim2. Calculates derivative at r, and takes step dr.
    */

    int nc = numc();
    int nq = numq();
    int i;

    double *dprim = (double *) malloc(nc * sizeof(double));

    flow_grad(prim1, r, dprim);
//    for(i=0; i<nc; i++)
//        printf("  %.12g", dprim[i]);
//    printf("\n");

    for(i=0; i<nc; i++)
        prim2[i] = prim1[i] + dr*dprim[i];
    for(i=nc; i<nq; i++)
        prim2[i] = prim1[i];

    free(dprim);
}

void substep_backward(double *prim1, double *prim2, double r, double dr)
{
    /*
    Primitive implicit first-order step.  Uses prim1 as input, saves result
    in prim2. Calculates derivative at r+dr iteratively, and takes step dr.
    */

    int max_count = 1000000;
    double tolerance = 1.0e-6;
    int nc = numc();
    int nq = numq();
    int i, j;

    double err = 0;
    double *dprim = (double *) malloc(nc * sizeof(double));
    double *old = (double *) malloc(nq * sizeof(double));
    double *new = (double *) malloc(nq * sizeof(double));

    flow_grad(prim1, r+dr, dprim);
    for(i=0; i<nc; i++)
    {    
        old[i] = prim1[i];
        new[i] = prim1[i] + dr*dprim[i];
        err += (new[i]-old[i]) * (new[i]-old[i]) / (old[i]*old[i]);
    }
    for(i=nc; i<nq; i++)
    {
        old[i] = prim1[i];
        new[i] = prim1[i];
    }

    j = 0;
    while( err > tolerance && j < max_count)
    {
        err = 0;
        for(i=0; i<nc; i++)
            old[i] = new[i];
        flow_grad(old, r+dr, dprim);
        for(i=0; i<nc; i++)
        {
            new[i] = prim1[i] + dr*dprim[i];
            err += (new[i]-old[i])*(new[i]-old[i]) / (old[i]*old[i]);
        }
        j++;
    }

    for(i=0; i<nq; i++)
        prim2[i] = new[i];

    printf("   Used %d iterations.\n", j);

    free(dprim);
    free(old);
    free(new);
}

void forward_euler(double *prim, double r, double *dr)
{
    /*
    First order explicit forward step.  
    */

    substep_forward(prim, prim, r, *dr);
}

void backward_euler(double *prim, double r, double *dr)
{
    /*
    First order implicit step.  
    */
    int i, over;
    int nq = numq();
    int nc = numc();
    double *prim2 = (double *) malloc(nq  * sizeof(double));
    do{
        substep_backward(prim, prim2, r, *dr);
        over = 0; 
        for(i=0; i<nc; i++)
        {
            if(fabs(prim2[i]) > 10.0*fabs(prim[i])
             || fabs(prim2[i]) < 0.1*fabs(prim[i]))
            {
                over = 1;
                *dr *= 0.5;
                printf("   dr: %.12lg prim:(%.12lg %.12lg %.12lg %.12lg)\n", 
                    *dr, prim[RHO], prim[TTT], prim[URR], prim[UPP]);
                break;
            }
        }
    }
    while(over);

    for(i=0; i<nq; i++)
        prim[i] = prim2[i];
    free(prim2);
}

void forward_rk2(double *prim, double r, double *dr)
{
    /*
    Second order explicit Runge-Kutta step.
    */

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
    /*
    Fourth order explicit Runge-Kutta step.
    */

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
