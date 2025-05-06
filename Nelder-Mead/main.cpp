#include "simpleks.hpp"

double funckja(double* punkt)
{
   return (pow((pow(punkt[0],2) - pow(punkt[1],2) - 1),2) + pow(punkt[1] - 1,2)+pow(punkt[2]-1,2));
   //return pow(punkt[0],2)+pow(punkt[1],2);//+pow(punkt[2],2);
   //return abs(punkt[0]+punkt[1]);//sin(punkt[1])+cos(pow(punkt[0],2))-punkt[0]*punkt[1];
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