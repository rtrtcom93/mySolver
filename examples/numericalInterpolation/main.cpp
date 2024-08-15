#include <iostream>
#include <cmath>
#include <fstream>
#include "interpolation.hpp"

using namespace std;
double fn(double x);
double *create_array(size_t n, double init_value=0);

int main(int argc, char** argv)
{
    const double PI{4*atan(1)};
    size_t nx{129};
    size_t num_data{11};
    double dx{0}, x0{-1}, xn{1};
    double *xdata{nullptr}, *zdata{nullptr}, *rdata{nullptr};
    double *ydata{nullptr}, *wdata{nullptr}, *sdata{nullptr}, *pdata{nullptr};

    dx = (xn - x0)/static_cast<double>(num_data-1);
    
    zdata = create_array(num_data);
    xdata = create_array(num_data);
    ydata = create_array(num_data);
    wdata = create_array(num_data); 
    rdata = create_array(nx);
    sdata = create_array(nx);
    pdata = create_array(nx);
    
    
    for (size_t i{0}; i<num_data; ++i)
    {
        xdata[i] = x0 + dx*static_cast<double>(i);
        zdata[i] = cos((2.*static_cast<double>(i)+1.)/(2.*static_cast<double>(num_data))*PI);
        ydata[i] = fn(xdata[i]);
        wdata[i] = fn(zdata[i]);
    }
    
    dx = (xn - x0)/static_cast<double>(nx-1);
    for (size_t i{0}; i < nx; ++i)
    {
        rdata[i] = x0 + dx*static_cast<double>(i);
        sdata[i] = lagrange(rdata[i], xdata, ydata, static_cast<int>(num_data));
        pdata[i] = lagrange(rdata[i], zdata, wdata, static_cast<int>(num_data));
    }

    ofstream Datafile("../result/data.txt");
    for (size_t i{0}; i < num_data; ++i)
    {
        Datafile << xdata[i] << "\t" << ydata[i] << "\t" << zdata[i] << "\t" << wdata[i] << endl;
    }
    Datafile.close();
   
    ofstream Resultfile("../result/result.txt");
    for (size_t i{0}; i < nx; ++i)
    {
        Resultfile << rdata[i] << "\t" << sdata[i] << "\t" << pdata[i] << endl;
    }
    Resultfile.close();

    delete [] zdata;
    delete [] xdata;
    delete [] ydata;
    delete [] wdata;
    delete [] sdata;
    delete [] rdata;
    delete [] pdata;

    return 0;
}

double fn(double x)
{
    return 1/(1+16*pow(x, 2));
}

double *create_array(size_t n, double init_value)
{
    double *new_array{nullptr};
    new_array = new double[n];
    for (size_t i{0}; i<n;++i)
        new_array[i] = init_value;
    return new_array;
}
