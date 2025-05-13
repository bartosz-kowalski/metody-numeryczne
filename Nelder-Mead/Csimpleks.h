#ifndef CSIMP
#define CSIMP
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../standard-calc/mean.h"
#include "../standard-calc/std.h"

typedef double (* funkcja)(double*);

void init(double**, const double*, const int);

double* calculate_mean_2d(double**, const int);

int compare(const void*, const void*);

void idx_sort(double*, double**, const int);

double* min_search(const int, double*, funkcja);
#endif
