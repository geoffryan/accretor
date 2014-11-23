#ifndef ACCRETOR_STEP
#define ACCRETOR_STEP

struct parList;

void (*step)(double *, double, double *);

void step_setup(struct parList *theParList);
void evolve(double *prim, double R1, double R2, struct parList *theParList);
void substep_forward(double *prim1, double *prim2, double r, double dr);
void substep_backward(double *prim1, double *prim2, double r, double dr);


void forward_euler(double *prim, double r, double *dr);
void backward_euler(double *prim, double r, double *dr);
void forward_rk2(double *prim, double r, double *dr);
void forward_rk4(double *prim, double r, double *dr);

#endif
