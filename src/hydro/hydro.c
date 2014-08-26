#include "hydro.h"

void hydro_setup(int choice)
{
    if(choice == 0)
    {
        numq = &numq_test;
        numc = &numc_test;
        initial = &initial_test;
        flow_grad = &flow_grad_test;
    }
    else if(choice == 1)
    {
        numq = &numq_rel_disc;
        numc = &numc_rel_disc;
        initial = &initial_rel_disc;
        flow_grad = &flow_grad_rel_disc;
    }
}
