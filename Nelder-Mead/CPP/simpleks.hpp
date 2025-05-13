#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <iostream>

#pragma once

typedef double (* funkcja)(double*);

class simplex;

double calculate_mean(const double* f_values, int dim);

double* calculate_mean_2d(double** simp, int f_dim);

double calculate_std(const double* f_values, int dim);

void idx_sort(double* f_values, simplex& simp);

void printSimp(simplex* simp, funkcja func);

double* min_search(const int no_dim, double* start_pt, funkcja function);
