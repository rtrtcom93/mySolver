#ifndef ROOT_FINDING_HPP
#define ROOT_FINDING_HPP

#include <iostream>
#include <cmath>

double bisect(double (*fn)(double), double x0, double xn, int maxIter=100, double tol=1e-12, bool prnt=true);
double bisect(double (*fn)(double x, int n), int order, double x0, double xn, int maxIter=100, double tol=1e-12, bool prnt=true);

double newton1D(double (*fn)(double), double(*dfn)(double), double init_guess, int maxIter=100, double tol=1e-12, bool prnt=true);
double secant1D(double (*fn)(double), double init_dx, double init_guess, double alpha=1., int maxIter=100, double tol=1e-12, bool prnt=true);

#endif