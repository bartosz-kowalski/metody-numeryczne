#include "gradient.h"

const double h = 0.000000001;

double* gradient(const int size, double* values, funkcja function)
{
    double* grad = (double*)malloc(sizeof(double)*size);
    double* f = (double*)malloc(sizeof(double)*size);
    double f1, f2, f3, f4;
    //roznica skonczona 4 rzedu

    for(int i = 0; i < size; ++i)
    {
        memcpy(f, values, size * sizeof(double));
        f[i] += 2*h;
        f1 = function(f);
        f[i] -= h;
        f2 = function(f);
        f[i] -= 2*h;
        f3 = function(f);
        f[i] -= h;
        f4 = function(f);
        grad[i]=(-f1+8*f2-8*f3+f4)/(12*h);
    }

    free(f);

    return grad;
}

