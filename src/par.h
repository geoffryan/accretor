#ifndef ACCRETOR_PAR
#define ACCRETOR_PAR

enum{VAR_DBL, VAR_INT, VAR_LON};

struct parList
{
    int hydro;
    int eos;
    int cool;
    int step;
    int N;
    int stop;
};

const static struct parList PAR_DEFAULT = {
    .hydro = 0,
    .eos = 0,
    .cool = 0,
    .step = 0,
    .N = 1,
    .stop = -1
};

int readvar(char filename[], char key[], int vtype, void *ptr);
void read_pars(struct parList *theParList, char filename[]);

#endif
