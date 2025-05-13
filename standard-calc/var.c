#include "var.h"

double variance(const double* values, const int size)
{
    double mean_ = mean(values, size);
    double variance = 0;
    for( int i = 0; i < size; ++i)
    {
        variance += pow((values[i] - mean_), 2);
    }
    return variance/(size-1);
} 
