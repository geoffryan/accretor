int numq_test()
{
    return 3;
}

int numc_test()
{
    return 2;
}

void initial_test(double *prim)
{
    prim[0] = 1.0;
    prim[1] = 0.0;
    prim[2] = 1.0;
}

void flow_grad_test(double *prim, double x, double *dprim)
{
    dprim[0] = prim[1];
    dprim[1] = prim[2]*prim[0];
}
