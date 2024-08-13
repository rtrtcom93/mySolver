#include <iostream>
#include "../include/root_finding.hpp"

double fn(double x);
double dfn(double x);

int main(int argc, char** argv)
{

    //Bisect(fn, -10, -1);
    //Bisect(fn, -1,   0);
    //Bisect(fn,  0,  10);

    newton1D(fn, dfn, -10);
    //Newton1D(fn, dfn, -0.1);
    //Newton1D(fn, dfn,  10);

    secant1D(fn, -0.01, -10);


    return 0; 
}

double fn(double x)
{
    return pow(x, 5) - 9*pow(x, 4) - pow(x, 3) + 17*pow(x, 2)- 8*x - 8;
}

double dfn(double x)
{
    return 5*pow(x, 4) - 36*pow(x, 3) - 3*pow(x, 2) + 34*x - 8;
}