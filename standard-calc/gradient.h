#ifndef GRADIENT
#define GRADIENT

#include <string.h>
#include <stdlib.h>

//const double h = 0.00000000001;

typedef double (* funkcja)(double*);

double* gradient(const int, double*, funkcja);

#endif
