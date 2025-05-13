#ifndef NEWTONC
#define NEWTONC

#include <stdlib.h>
#include <math.h>
#include "../standard-calc/gradient.h"

typedef double (* funkcja)(double*);

double* newton_raphson(const int, double*, funkcja);

#endif