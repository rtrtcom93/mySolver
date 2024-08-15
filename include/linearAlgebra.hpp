#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

void displayVector(const vector<double> &vec);
void displayMatrix(const vector<vector<double>> &matrix);

double norm(const vector<double> &vec, string &&type);
double dotProd(const vector<double> &vec1, const vector<double> &vec2);

void replace(vector<double> &vec1, vector<double> &vec2, size_t m);
void replace(vector<double> &vec1, vector<double> &&vec2, size_t m);
vector<double> slice(const vector<double> &vec, int m, int n);
vector<double> linspace(double x0, double xn, size_t size);

vector<double> matMul(const vector<vector<double>>& mat, const vector<double> vec);

void gaussJordan(vector<vector<double>>& matrix);
void gaussJordanInv(vector<vector<double>>& matrix);
vector<double> TDMA(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d);

#endif