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
}
