#include "../include/rootFinding.hpp"

double bisect(double (*fn)(double), double x0, double xn, int maxIter, double tol, bool prnt)
{
    using namespace std;

    double res{1.};
    double xm{0.};

    cout.precision(12);
    if (prnt == true)
    {
        cout << "iteration\t" << "residual\t" << "x0\t" << "xm\t" << "xn\n";
    }
    for (int i{0}; i<maxIter; i++)
    {   
        xm = (xn + x0)/2.;
        res = abs(xn - xm);
        if (prnt == true)
            cout << i << "\t"<< res << "\t" << x0 << "\t" << xm << "\t" << xn << endl;
        if (res < tol) 
            break;
        else
            (fn(xm)*fn(xn) <= 0) ? x0 = xm : xn = xm;  
    }

    return xm;
}

double bisect(double (*fn)(double x, int n), int order, double x0, double xn, int maxIter, double tol, bool prnt)
{
    using namespace std;

    double res{1.};
    double xm{0.};

    cout.precision(12);
    if (prnt == true)
    {
        cout << "iteration\t" << "residual\t" << "x0\t" << "xm\t" << "xn\n";
    }
    for (int i{0}; i<maxIter; i++)
    {   
        xm = (xn + x0)/2.;
        res = abs(xn - xm);
        if (prnt == true)
            cout << i << "\t"<< res << "\t" << x0 << "\t" << xm << "\t" << xn << endl;
        if (res < tol) 
            break;
        else
            (fn(xm, order)*fn(xn, order) <= 0) ? x0 = xm : xn = xm;  
    }

    return xm;
}

double newton1D(double (*fn)(double), double(*dfn)(double), double init_guess, int maxIter, double tol, bool prnt)
{
    using namespace std;
    
    double xb{init_guess}, xn{0.};
    double res{1.};

    cout.precision(12);
    if (prnt == true)
    {
        cout << "iteration\t" << "residual\t" << "fn\t" << "dfn\t" << "x\n" ;
    }

    xb = init_guess;
    for (int i{0}; i<maxIter; i++)
    {   
        if (dfn(xb) == 0) 
            cout << "Zero division" << endl;
        xn = xb-fn(xb)/dfn(xb);
        res = abs(xn - xb);
        if (prnt == true)
            cout << i << "\t"<< res << "\t" << fn(xb) << "\t" << dfn(xb) << "\t" << xb << endl;     
        if (res < tol)  
            break;
        else 
            xb = xn;
    }

    return xn;
}

double secant1D(double (*fn)(double), double init_dx, double init_guess, double alpha, int maxIter, double tol, bool prnt)
{
    using namespace std;
    
    double dfn{0};
    double xb{init_guess}, xn{0.};

    cout.precision(12);
    if (prnt == true)
    {
        cout << "iteration\t" << "residual\t" << "fn\t" << "dfn\t" << "x\n" ;
    }

    xb = init_guess;
    dfn = (fn(init_guess+init_dx)-fn(init_guess))/init_dx;
    for (int i{0}; i < maxIter; i++)
    {   
        xn = xb - alpha*fn(xb)/dfn;
        if (prnt == true)
            cout << i << "\t"<< abs(xn - xb) << "\t" << fn(xb) << "\t" << dfn << "\t" << xb << endl;     
        if (abs(xn - xb)< tol)  
            break;
        else if (abs(fn(xn)) < tol)
            break;
        else {
            dfn = (fn(xn) - fn(xb))/(xn - xb);
            xb = xn;
        }
    }

    return xn;
}