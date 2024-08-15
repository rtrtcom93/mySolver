#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include "linearAlgebra.hpp"
using namespace std;

double lagrange_i(double x, unsigned i, double *xdata, unsigned numData);
double lagrange_i(double x, unsigned i, const vector<double> &xdata);
double lagrange(double x, const double *const xdata, const double *const ydata, int numData);
double lagrange(double x, const vector<double> &xdata, const vector<double> &ydata);

double cubicSpline(double x, const vector<double>& xdata, const vector<double>& ydata, string bc_type="natural", double para_a=1., double para_b=1.);
vector<double> cubicSpline(vector<double> x, const vector<double>& xdata, const vector<double>& ydata, string bc_type="natural", double para_a=1., double para_b=1.);

#endif