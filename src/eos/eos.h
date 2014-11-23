#ifndef ACCRETOR_EOS
#define ACCRETOR_EOS

struct parList;

static double GAMMA = 4.0/3.0;
static double visc_f = 1.0;
static double q0 = 1.0e22;

double (*pressure)(double *, double);
double (*spec_int_en)(double *, double);
double (*depsdrho)(double *, double);
double (*depsdT)(double *, double);
double (*dPdrho)(double *, double);
double (*dPdT)(double *, double);
double (*cool)(double *, double);

void eos_setup(struct parList *theParList);
void cool_setup(struct parList *theParList);

double pressure_ideal(double *prim, double r);
double spec_int_en_ideal(double *prim, double r);
double depsdrho_ideal(double *prim, double r);
double depsdT_ideal(double *prim, double r);
double dPdrho_ideal(double *prim, double r);
double dPdT_ideal(double *prim, double r);

double cool_none(double *prim, double r);
double cool_thompson(double *prim, double r);
double cool_visc(double *prim, double r);


#endif
