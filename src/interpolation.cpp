#include "../include/interpolation.hpp"

double lagrange_i(double x, unsigned i, const double *const xdata, unsigned numData)
{
    double val{1.};
    
    for (unsigned j{0}; j<numData; ++j) 
        (j == i) ? val : val *= (x - xdata[j])/(xdata[i] - xdata[j]) ;        

    return val;
}

double lagrange_i(double x, unsigned i, const vector<double> &xdata)
{
    size_t numData{xdata.size()};
    double val{1.};
    
    for (unsigned j{0}; j < numData; ++j) 
        (j == i) ? val : val *= (x - xdata[j])/(xdata[i] - xdata[j]) ;        

    return val;
}

double lagrange(double x, const double *const xdata, const double *const ydata, int numData)
{
    double val{0.};
    
    for (int i{0}; i < numData; ++i)
        val += ydata[i]*lagrange_i(x, i, xdata, numData);
    
    return val; 
}

double lagrange(double x, const vector<double> &xdata, const vector<double> &ydata)
{
    size_t numData{xdata.size()};
    double val{0.};
    
    for (unsigned i{0}; i < numData; ++i)
        val += ydata[i]*lagrange_i(x, i, xdata);
    
    return val; 
}

double cubicSpline(double x, const vector<double>& xdata, const vector<double>& ydata, string bc_type, double para_a, double para_b)
{
    double g_ix{0};

    size_t n1{xdata.size()};
    size_t n1m{n1-1};
    
    vector<double> dx(n1);
    vector<double> a(n1), b(n1), c(n1), rhs(n1);
    vector<double> ddg(n1);


    enum BoundaryType {
        natural, parabolic, periodic
    };

    BoundaryType bc_type_enum;

    if (bc_type == "parabolic") {
        bc_type_enum = parabolic;
    } else if (bc_type == "periodic")
    {
        bc_type_enum = periodic;
    } else {
        bc_type_enum = natural;
    }
    
    for (size_t i{0}; i < n1m; ++i)
        dx.at(i) = xdata.at(i+1)-xdata.at(i);     

 

    //Seting up for Tri-diagoanl matrix
    switch (bc_type_enum) {
    case parabolic :
        c.at(1)   = dx.at(1)/6.;
        b.at(1)   = (dx.at(1) + dx.at(0))/3.
                  + para_a*dx.at(0)/6.;
        rhs.at(1) = (ydata.at(2) - ydata.at(1))/dx.at(1) 
                - (ydata.at(1) - ydata.at(0))/dx.at(0);
        
        for (size_t i{2}; i < n1m-1; ++i)
            a.at(i) = dx.at(i-1)/6.;
        for (size_t i{2}; i < n1m-1; ++i) {
            b.at(i) = (dx.at(i)+dx.at(i-1))/3.;
            rhs.at(i) = (ydata.at(i+1) - ydata.at(i))/dx.at(i) 
                      - (ydata.at(i) - ydata.at(i-1))/dx.at(i-1);
        }
        for (size_t i{2}; i <n1m-1;++i)
            c.at(i) = dx.at(i+1)/6.;  
        
        a.at(n1m-1)   = dx.at(n1m-2)/6.;
        b.at(n1m-1)   = (dx.at(n1m-1) + dx.at(n1m-2))/3.
                      + para_b*dx.at(n1m-1)/6.;
        rhs.at(n1m-1) = (ydata.at(n1m) - ydata.at(n1m-1))/dx.at(n1m-1) 
                    - (ydata.at(n1m-1) - ydata.at(n1m-2))/dx.at(n1m-2);
        
        //solving TDMA
        rhs = TDMA(slice(a, 1, n1m), slice(b, 1, n1m), slice(c, 1, n1m), slice(rhs, 1, n1m));
        
        for (size_t i{1}; i < n1m; ++i)
            ddg[i] = rhs[i];  

        break;

    case periodic :
        c.at(1)   = dx.at(1)/6.;
        b.at(1)   = (dx.at(1) + dx.at(0))/3.;
        rhs.at(1) = (ydata.at(2) - ydata.at(1))/dx.at(1) 
                - (ydata.at(1) - ydata.at(0))/dx.at(0);
        
        for (size_t i{2}; i < n1m-1; ++i)
            a.at(i) = dx.at(i-1)/6.;
    
        for (size_t i{2}; i < n1m-1; ++i) {
            b.at(i) = (dx.at(i)+dx.at(i-1))/3.;
            rhs.at(i) = (ydata.at(i+1) - ydata.at(i))/dx.at(i) 
                      - (ydata.at(i) - ydata.at(i-1))/dx.at(i-1);
        }
        for (size_t i{2}; i <n1m-1;++i)
            c.at(i) = dx.at(i+1)/6.;

        a.at(n1m-1)   = dx.at(n1m-2)/6.;
        b.at(n1m-1)   = (dx.at(n1m-1) + dx.at(n1m-2))/3.
                      + para_b*dx.at(n1m-1)/6.;
        rhs.at(n1m-1) = (ydata.at(n1m) - ydata.at(n1m-1))/dx.at(n1m-1) 
                    - (ydata.at(n1m-1) - ydata.at(n1m-2))/dx.at(n1m-2);

        //solving TDMA
        rhs = TDMA(slice(a, 1, n1m), slice(b, 1, n1m), slice(c, 1, n1m), slice(rhs, 1, n1m));
        
        for (size_t i{1}; i < n1m; ++i)
            ddg[i] = rhs[i];

        ddg[0]   = ddg[n1m-1];
        ddg[n1m] = ddg[1];  
        
        break;

    default:
        c.at(1)   = dx.at(1)/6.;
        b.at(1)   = (dx.at(1) + dx.at(0))/3.;
        rhs.at(1) = (ydata.at(2) - ydata.at(1))/dx.at(1) 
                - (ydata.at(1) - ydata.at(0))/dx.at(0);  
        
        for (size_t i{2}; i < n1m-1; ++i)
            a.at(i) = dx.at(i-1)/6.;
    
        for (size_t i{2}; i < n1m-1; ++i) {
            b.at(i) = (dx.at(i)+dx.at(i-1))/3.;
            rhs.at(i) = (ydata.at(i+1) - ydata.at(i))/dx.at(i) 
                      - (ydata.at(i) - ydata.at(i-1))/dx.at(i-1);
        }
        for (size_t i{2}; i <n1m-1;++i)
            c.at(i) = dx.at(i+1)/6.;
        
        a.at(n1m-1)   = dx.at(n1m-2)/6.;
        b.at(n1m-1)   = (dx.at(n1m-1) + dx.at(n1m-2))/3.;
        rhs.at(n1m-1) = (ydata.at(n1m) - ydata.at(n1m-1))/dx.at(n1m-1) 
                    - (ydata.at(n1m-1) - ydata.at(n1m-2))/dx.at(n1m-2);

        //solving TDMA
        rhs = TDMA(slice(a, 1, n1m), slice(b, 1, n1m), slice(c, 1, n1m), slice(rhs, 1, n1m));  

        for (size_t i{1}; i < n1m; ++i)
            ddg[i] = rhs[i];

        break;
    };

//    displayVector(ddg);

    size_t r{0};
    if (x > xdata[n1m] || x < xdata[0]) {
        cout << "The point is out-of-range\n";
    } else {
        while (x > xdata[r]) {
            ++r;
        }
    }
    g_ix = ddg[r]/6.*(pow((xdata[r  ] - x), 3)/dx[r-1] - dx[r-1]*(xdata[r  ] - x))
         + ddg[r]/6.*(pow((x - xdata[r-1]), 3)/dx[r-1] - dx[r-1]*(x - xdata[r-1]))
         + ydata[r-1]*(xdata[r] - x)/dx[r-1]
         + ydata[r]*(x - xdata[r-1])/dx[r-1];

    return g_ix;
}

vector<double> cubicSpline(vector<double> x, const vector<double>& xdata, const vector<double>& ydata, string bc_type, double para_a, double para_b)
{
    vector<double> g_ix(x.size());

    size_t n1{xdata.size()};
    size_t n1m{n1-1};
    
    vector<double> dx(n1);
    vector<double> a(n1), b(n1), c(n1), rhs(n1);
    vector<double> ddg(n1);


    enum BoundaryType {
        natural, parabolic, periodic
    };

    BoundaryType bc_type_enum;

    if (bc_type == "parabolic") {
        bc_type_enum = parabolic;
    } else if (bc_type == "periodic")
    {
        bc_type_enum = periodic;
    } else {
        bc_type_enum = natural;
    }
    
    for (size_t i{0}; i < n1m; ++i)
        dx.at(i) = xdata.at(i+1)-xdata.at(i);     

 

    //Seting up for Tri-diagoanl matrix
    switch (bc_type_enum) {
    case parabolic :
        c.at(1)   = dx.at(1)/6.;
        b.at(1)   = (dx.at(1) + dx.at(0))/3.
                  + para_a*dx.at(0)/6.;
        rhs.at(1) = (ydata.at(2) - ydata.at(1))/dx.at(1) 
                - (ydata.at(1) - ydata.at(0))/dx.at(0);
        
        for (size_t i{2}; i < n1m-1; ++i)
            a.at(i) = dx.at(i-1)/6.;
        for (size_t i{2}; i < n1m-1; ++i) {
            b.at(i) = (dx.at(i)+dx.at(i-1))/3.;
            rhs.at(i) = (ydata.at(i+1) - ydata.at(i))/dx.at(i) 
                      - (ydata.at(i) - ydata.at(i-1))/dx.at(i-1);
        }
        for (size_t i{2}; i <n1m-1;++i)
            c.at(i) = dx.at(i+1)/6.;  
        
        a.at(n1m-1)   = dx.at(n1m-2)/6.;
        b.at(n1m-1)   = (dx.at(n1m-1) + dx.at(n1m-2))/3.
                      + para_b*dx.at(n1m-1)/6.;
        rhs.at(n1m-1) = (ydata.at(n1m) - ydata.at(n1m-1))/dx.at(n1m-1) 
                    - (ydata.at(n1m-1) - ydata.at(n1m-2))/dx.at(n1m-2);
        
        //solving TDMA
        rhs = TDMA(slice(a, 1, n1m), slice(b, 1, n1m), slice(c, 1, n1m), slice(rhs, 1, n1m));
        
        for (size_t i{1}; i < n1m; ++i)
            ddg[i] = rhs[i];  

        break;

    case periodic :
        c.at(1)   = dx.at(1)/6.;
        b.at(1)   = (dx.at(1) + dx.at(0))/3.;
        rhs.at(1) = (ydata.at(2) - ydata.at(1))/dx.at(1) 
                - (ydata.at(1) - ydata.at(0))/dx.at(0);
        
        for (size_t i{2}; i < n1m-1; ++i)
            a.at(i) = dx.at(i-1)/6.;
    
        for (size_t i{2}; i < n1m-1; ++i) {
            b.at(i) = (dx.at(i)+dx.at(i-1))/3.;
            rhs.at(i) = (ydata.at(i+1) - ydata.at(i))/dx.at(i) 
                      - (ydata.at(i) - ydata.at(i-1))/dx.at(i-1);
        }
        for (size_t i{2}; i <n1m-1;++i)
            c.at(i) = dx.at(i+1)/6.;

        a.at(n1m-1)   = dx.at(n1m-2)/6.;
        b.at(n1m-1)   = (dx.at(n1m-1) + dx.at(n1m-2))/3.
                      + para_b*dx.at(n1m-1)/6.;
        rhs.at(n1m-1) = (ydata.at(n1m) - ydata.at(n1m-1))/dx.at(n1m-1) 
                    - (ydata.at(n1m-1) - ydata.at(n1m-2))/dx.at(n1m-2);

        //solving TDMA
        rhs = TDMA(slice(a, 1, n1m), slice(b, 1, n1m), slice(c, 1, n1m), slice(rhs, 1, n1m));
        
        for (size_t i{1}; i < n1m; ++i)
            ddg[i] = rhs[i];

        ddg[0]   = ddg[n1m-1];
        ddg[n1m] = ddg[1];  
        
        break;

    default:
        c.at(1)   = dx.at(1)/6.;
        b.at(1)   = (dx.at(1) + dx.at(0))/3.;
        rhs.at(1) = (ydata.at(2) - ydata.at(1))/dx.at(1) 
                - (ydata.at(1) - ydata.at(0))/dx.at(0);  
        
        for (size_t i{2}; i < n1m-1; ++i)
            a.at(i) = dx.at(i-1)/6.;
    
        for (size_t i{2}; i < n1m-1; ++i) {
            b.at(i) = (dx.at(i)+dx.at(i-1))/3.;
            rhs.at(i) = (ydata.at(i+1) - ydata.at(i))/dx.at(i) 
                      - (ydata.at(i) - ydata.at(i-1))/dx.at(i-1);
        }
        for (size_t i{2}; i <n1m-1;++i)
            c.at(i) = dx.at(i+1)/6.;
        
        a.at(n1m-1)   = dx.at(n1m-2)/6.;
        b.at(n1m-1)   = (dx.at(n1m-1) + dx.at(n1m-2))/3.;
        rhs.at(n1m-1) = (ydata.at(n1m) - ydata.at(n1m-1))/dx.at(n1m-1) 
                    - (ydata.at(n1m-1) - ydata.at(n1m-2))/dx.at(n1m-2);

        //solving TDMA
        rhs = TDMA(slice(a, 1, n1m), slice(b, 1, n1m), slice(c, 1, n1m), slice(rhs, 1, n1m));  

        for (size_t i{1}; i < n1m; ++i)
            ddg[i] = rhs[i];

        break;
    };

//    displayVector(ddg);

    size_t r{0};
    for (size_t i{0}; i < x.size(); ++i) {
        r = 0;
        if (x[i] > xdata[n1m] || x[i] < xdata[0]) {
            cout << "The point is out-of-range\n";
        } else {
            while (x[i] > xdata[r]) {
                ++r;
            }
        }
        g_ix[i] = ddg[r]/6.*(pow((xdata[r] - x[i]), 3)/dx[r-1] - dx[r-1]*(xdata[r]-x[i]))
                + ddg[r]/6.*(pow((x[i]-xdata[r-1]), 3)/dx[r-1] - dx[r-1]*(x[i]-xdata[r-1]))
                + ydata[r-1]*(xdata[r] - x[i])/dx[r-1]
                + ydata[r]*(x[i] - xdata[r-1])/dx[r-1];
    }

    return g_ix;
}