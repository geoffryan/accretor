#ifndef ACCRETOR_HYDRO
#define ACCRETOR_HYDRO

#define RHO 0
#define TTT 1
#define URR 2
#define UPP 3
#define DUR 4
#define DUP 5
#define ACC 4
#define LLL 5
#define SAC 3

static double M = 1.0;
static double alpha = 0.01;

int (*numq)();
int (*numc)();
void (*initial)(double *, double *, double *);
void (*exact)(double *, double);
void (*flow_grad)(double *, double, double *);

void hydro_setup(int choice);

int numq_test();
int numc_test();
void initial_test(double *prim, double *R1, double *R2);
void exact_test(double *prim, double x);
void flow_grad_test(double *prim, double x, double *dprim);
int numq_rel_disc();
int numc_rel_disc();
void initial_rel_disc(double *prim, double *R1, double *R2);
void exact_rel_disc(double *prim, double r);
void flow_grad_rel_disc(double *prim, double x, double *dprim);
int numq_newt_disc_SS();
int numc_newt_disc_SS();
void initial_newt_disc_SS(double *prim, double *R1, double *R2);
void exact_newt_disc_SS(double *prim, double r);
void flow_grad_newt_disc_SS(double *prim, double x, double *dprim);
int numq_newt_sph();
int numc_newt_sph();
void initial_newt_sph(double *prim, double *R1, double *R2);
void exact_newt_sph(double *prim, double r);
void flow_grad_newt_sph(double *prim, double x, double *dprim);

#endif
