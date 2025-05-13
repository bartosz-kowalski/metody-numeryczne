#ifndef CSIMPM
#define CSIMPM

#include "Csimpleks.h"

#define ALPHA 1.1
#define BETA 0.4
#define GAMMA 2.05
#define DELTA 0.3
#define EPSILON 0.00001
#define MAX_I 1000
#define SIZE(X) sizeof(X)/sizeof(typeof(X[0]))

void init(double** simplex, const double* starting_point, const int dim)
{   
    //inicjalizacja simpleksu
    for(int i = 0; i < dim; ++i)
    {
        simplex[0][i]=starting_point[i];
    }

    //uzupelnienie o inne punkty
    for(int j = 1; j < dim+1; ++j)
    {
        for(int i = 0, k = 1; i < dim; ++i, ++k)
        {
            simplex[j][i]=simplex[0][i];
            simplex[k][i]=simplex[0][i]+0.1;
        }
    }
    return;
}

double* calculate_mean_2d(double** simp, const int f_dim)
{
    double* mean =  (double*)malloc(sizeof(double)*f_dim);

    for(int i = 0; i < f_dim; ++i)
    {   
        mean[i] = 0;
        for(int j = 0; j < f_dim; ++j)
        {
            mean[i] += simp[j][i];
        }
        mean[i] /= f_dim;
    }
    return mean;
}

int compare(const void *a, const void *b) {
    double i = *(const double *)a;
    double j = *(const double *)b;
    return (i > j) - (i < j);
}

void idx_sort(double* f_values, double** sim, const int dimensions)
{
    int dim_s = dimensions+1;    // liczba punktów w simpleksie
    int dim = dimensions;       // wymiar przestrzeni

    double* f_unsorted = (double*)malloc(sizeof(double)*dim_s);
    int* idx = (int*)malloc(sizeof(int)*dim_s);
    for (int i = 0; i < dim_s; ++i)
    {
        f_unsorted[i]=f_values[i];
    }
    // Posortuj indeksy względem f_values
    qsort(f_values, dim_s, sizeof(double), compare);

    for(int i = 0; i < dim_s; ++i)
    {
        for( int j = 0; j < dim_s; ++j)
        {
            if(f_values[i]==f_unsorted[j])
            {
                idx[i]=j;
            }
        }
    }

    double** sim_sorted = (double**)malloc(sizeof(double)*dim_s);
    for (int i = 0; i < dim_s; ++i) {
        sim_sorted[i] = (double*)malloc(sizeof(double)*dim);
        for (int j = 0; j < dim; ++j) {
            sim_sorted[i][j] = sim[idx[i]][j];
        }
    }

    for (int i = 0; i < dim_s; ++i) {
        for (int j = 0; j < dim; ++j) {
            sim[i][j] = sim_sorted[i][j];
        }
        free(sim_sorted[i]);
    }
    free(sim_sorted);
    free(idx);
    free(f_unsorted);
}

double* min_search(const int no_dim, double* start_pt, funkcja function)
{
    //int no_dim = sizeof(*start_pt)/sizeof(double);
    double** simplex;
    int simplex_point_no = no_dim+1;
    int dimension_no = no_dim;
   
    simplex = (double**)malloc(sizeof(double*)*simplex_point_no);
    for(int i = 0; i < simplex_point_no; ++i)
    {
        simplex[i] = (double*)malloc(sizeof(double)*dimension_no);
    }

    init(simplex, start_pt, dimension_no);

    //inicjalizacja wektora funkcji
    double* f_values = (double*)malloc(sizeof(double)*simplex_point_no);//new double[starting_point.getDimS()];
    for(int i = 0; i < simplex_point_no; ++i)
    {
        f_values[i]=function(simplex[i]);
    }

    idx_sort(f_values, simplex, dimension_no);
    
    //zmienne pomocnicze
    double* x_o = (double*)malloc(sizeof(double)*dimension_no);
    double* x_r = (double*)malloc(sizeof(double)*dimension_no);
    double* x_c = (double*)malloc(sizeof(double)*dimension_no);
    double* x_e = (double*)malloc(sizeof(double)*dimension_no);
    double f_r, f_c, f_e;

    //glowna pętla
    for(int i = 0; i < MAX_I; ++i)
    {
        //warunek wyjscia (nie dziala)
        if(std(f_values, simplex_point_no) < EPSILON) 
        {
            //printf("%s", "osiagnieto pozadane odchylenie standradowe \n");
            break;
        }

        idx_sort(f_values, simplex, dimension_no);

        x_o = calculate_mean_2d(simplex, dimension_no);

        for(int j = 0; j < dimension_no; ++j)
        {
            x_r[j] = x_o[j] + ALPHA * (x_o[j] - simplex[simplex_point_no-1][j]);
        }
        f_r=function(x_r);

        if(f_r >= f_values[0] && f_r < f_values[simplex_point_no-1])
        {
            for (int j = 0; j < dimension_no; ++j)
                simplex[simplex_point_no-1][j] = x_r[j];
            f_values[simplex_point_no-1] = f_r;
            continue;
        }
        else if(f_r < f_values[0])
        {
            for(int j = 0; j < dimension_no; ++j)
            {
                x_e[j] = x_o[j] + GAMMA * (x_r[j] - x_o[j]);
            }
            f_e=function(x_e);
            
            if(f_e < f_r)
            {
                for (int j = 0; j < dimension_no; ++j)
                    simplex[simplex_point_no-1][j] = x_e[j];
                f_values[simplex_point_no-1] = f_e;
                continue;
            }
            else
            {
                for (int j = 0; j < dimension_no; ++j)
                    simplex[simplex_point_no-1][j] = x_r[j];
                f_values[simplex_point_no-1] = f_r;
                continue;
            }
        }
        else if(f_r >= f_values[simplex_point_no-1])
        {
            for(int j = 0; j < dimension_no; ++j)
            {
                x_c[j] = x_o[j] + BETA * (simplex[simplex_point_no-1][j] - x_o[j]);
            }
            f_c=function(x_c);
            
            if(f_c < f_r)
            {
                for (int j = 0; j < dimension_no; ++j)
                    simplex[simplex_point_no-1][j] = x_c[j];
                f_values[simplex_point_no-1] = f_c;
                continue;
            }
        }

        for(int j = 1; j < simplex_point_no; ++j)
        {   
            for(int k = 0; k < dimension_no; ++k)
            {
                simplex[j][k] = simplex[0][k] + DELTA * (simplex[j][k] - simplex[0][k]);
            }
            f_values[j] = function(simplex[j]);
        }
    }

    free(x_o);
    free(x_r);
    free(x_c);
    free(x_e);
    free(f_values);
    
    for(int i = simplex_point_no-1; i > 0; --i)
        free(simplex[i]);
    
    return simplex[0];
}

#endif