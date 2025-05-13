#include "mean.h"

double mean(const double* values, const int size)
{
    double mean_ = 0;
    for(int i = 0; i < size; ++i)
    {
        mean_ += values[i];
    }
    return mean_/size;
}
