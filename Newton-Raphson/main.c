#include "newton.h"
#include <stdio.h>

double funckja(double* punkt)
{
   //return (100*pow(-pow(punkt[0],2)+punkt[1],2)+pow((1-punkt[0]),2));
   return pow(punkt[0]-2,2)-30;
}

int main()
{
    int cols = 1;
    double* data = (double*)malloc(sizeof(double)*cols);

    data[0]=-1.5;
    //data[1]=1;
    //data[2]=-0.5;

    double* zerowe = newton_raphson(cols,data,funckja);
    printf("Miejsce zerowe funkcji: %6.10f, \n",zerowe[0]);//,minimum[1]);
    printf("Wartosc funkcji miejscu zerowym (xD): %6.8f \n",funckja(zerowe));

    //free(data);
    free(zerowe);

    return 0;
}