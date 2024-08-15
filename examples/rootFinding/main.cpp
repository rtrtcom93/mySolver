#include <iostream>
#include "rootFinding.hpp"

double fn(double x);
double dfn(double x);

int main(int argc, char** argv)
{
    using namespace std;

    cout << "Results : Newton method" << endl;
    newton1D(fn, dfn, -10.0);

    cout << "Results : Secant method" << endl;
    secant1D(fn, -0.01, -10.0);


    return 0; 
}

/*Gvien function f(x)
  f(x) = x^5 - 9x^4 - x^3 + 17x^2 - 8x - 8
  df(x)/dx = 5x^4 - 36x^3 - 3x^2 + 18*/
double fn(double x)
{
    return pow(x, 5) - 9*pow(x, 4) - pow(x, 3) + 17*pow(x, 2)- 8*x - 8;
}

double dfn(double x)
{
    return 5*pow(x, 4) - 36*pow(x, 3) - 3*pow(x, 2) + 34*x - 8;
}