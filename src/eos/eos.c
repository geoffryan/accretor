#include "../par.h"
#include "eos.h"

void eos_setup(struct parList *theParList)
{
    if(theParList->eos == 0)
    {
        pressure = &pressure_ideal;
        spec_int_en = &spec_int_en_ideal;
        dPdrho = &dPdrho_ideal;
        dPdT = &dPdT_ideal;
        depsdrho = &depsdrho_ideal;
        depsdT = &depsdT_ideal;
    }
}

void cool_setup(struct parList *theParList)
{
    if(theParList->cool == 0)
    {
        cool = &cool_none;
    }
    else if(theParList->cool == 1)
    {
        cool = &cool_thompson;
    }
    else if(theParList->cool == 2)
    {
        cool = &cool_visc;
    }
}
