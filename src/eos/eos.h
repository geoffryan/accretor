#ifndef ACCRETOR_EOS
#define ACCRETOR_EOS

double (*pressure)(double *, double);
double (*spec_int_en)(double *, double);
double (*depsdrho)(double *, double);
double (*depsdT)(double *, double);
double (*dPdrho)(double *, double);
double (*dPdT)(double *, double);
double (*cool)(double *, double);

void eos_setup(int choice);
void cool_setup(int choice);

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
