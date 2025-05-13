#include "gradient.h"
#include <stdio.h>
#include <math.h>

double funckja(double* punkt)
{
   return (100*pow(-pow(punkt[0],2)+punkt[1],2)+pow((1-punkt[0]),2));
   //return punkt[0];
}

int main()
{
    int cols = 2;
    double* data = (double*)malloc(sizeof(double)*cols);

    data[0]=-1.5;
    data[1]=1;
    //data[2]=-0.5;

    double* zerowe = gradient(cols,data,funckja);
    printf("Miejsce zerowe funkcji: %6.10f, %6.10f \n",zerowe[0], zerowe[1]);

    //free(data);
    free(zerowe);

    return 0;
}
