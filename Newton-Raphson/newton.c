#include "newton.h"

double* newton_raphson(const int dim, double* starting_point, funkcja function)
{
    const double eps = 0.0001;
    static int i = 0;

    double value = function(starting_point);

    if(i>5000)
        return starting_point;

    if(value < eps && value > -eps) 
        return starting_point;

    double* next_point = (double*)malloc(sizeof(double)*dim);

    double* grad = gradient(dim, starting_point, function);

    for(int i = 0; i < dim; ++i)
    {
        next_point[i] = starting_point[i] - value/grad[i];
    }

    free(grad);
    i++;
    return newton_raphson(dim, next_point, function);
}
