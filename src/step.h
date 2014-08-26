#ifndef ACCRETOR_STEP
#define ACCRETOR_STEP

void (*step)(double *, double, double *);

void step_setup(int choice);
void evolve(double *prim, double R1, double R2, int n);
void substep_forward(double *prim1, double *prim2, double r, double dr);

void forward_euler(double *prim, double r, double *dr);
void forward_rk2(double *prim, double r, double *dr);
void forward_rk4(double *prim, double r, double *dr);

#endif
