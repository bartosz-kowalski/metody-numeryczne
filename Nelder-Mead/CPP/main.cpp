#include "simpleks.hpp"

double funckja(double* punkt)
{
    return pow(punkt[0],2)+pow(punkt[1],2)+pow(punkt[2],2);
   //return (pow((pow(punkt[0],2) - pow(punkt[1],2) - 1),2) + pow(punkt[1] - 1,2)+pow(punkt[2]-1,2));
   //return (100*pow(-pow(punkt[0],2)+punkt[1],2)+pow((1-punkt[0]),2));
   //return (pow(1.5-punkt[0]+punkt[0]*punkt[1],2)+pow(2.25-punkt[0]-pow(punkt[0]*punkt[1],2),2)+pow(2.625-punkt[0]+pow(punkt[0]*punkt[1],3),2));
}

int main()
{
    int cols = 3;
    double* data = new double[cols];

    data[0]=2;
    data[1]=1;
    data[2]=-0.5;

    double* minimum = min_search(3,data,funckja);
    printf("Minimum funkcji: %6.5f, %6.5f, %6.5f \n",minimum[0],minimum[1], minimum[2]);
    printf("Wartosc funkcji w minimum: %6.8f \n",funckja(minimum));

    delete[] data;
    delete[] minimum;

    return 0;
}