#include "eos.h"

void eos_setup(int choice)
{
    if(choice == 0)
    {
        pressure = &pressure_ideal;
        spec_int_en = &spec_int_en_ideal;
        dPdrho = &dPdrho_ideal;
        dPdT = &dPdT_ideal;
        depsdrho = &depsdrho_ideal;
        depsdT = &depsdT_ideal;
    }
}

void cool_setup(int choice)
{
    if(choice == 0)
    {
        cool = &cool_none;
    }
    if(choice == 1)
    {
        cool = &cool_thompson;
    }
}
