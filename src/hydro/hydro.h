#ifndef ACCRETOR_HYDRO
#define ACCRETOR_HYDRO

int (*numq)();
int (*numc)();
void (*initial)(double *);
void (*flow_grad)(double *, double, double *);

void hydro_setup(int choice);

int numq_test();
int numc_test();
void initial_test(double *prim);
void flow_grad_test(double *prim, double x, double *dprim);

#endif
