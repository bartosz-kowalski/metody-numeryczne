#include "simpleks.hpp"

#define ALPHA 1.1
#define BETA 0.4
#define GAMMA 2.05
#define DELTA 0.3
#define EPSILON 0.0000001
#define MAX_I 1000


class simplex{
    private:
        const int dim; //number of dimensions
        const int dim_sim;
        double** sim; // the simplex itself

    public:
        simplex();
        simplex(int);
        ~simplex();

        void update(double**); //updates the simplex's dimensions and the simplex coordinates
        double** getSimp(); //returns the simplex
        int getDim(); //returns the no. of dimensions
        int getDimS(); //returns the no. of dimensions of the simplex (one more then getDim())
        void init(double*); //initializes the simplex
    
};

simplex::simplex()
    : dim(0)
    ,dim_sim(0)
    , sim(nullptr)
{
}

simplex::simplex(int num)
    : dim(num)
    , dim_sim(num+1)
{
    sim = new double*[dim_sim];

    for(int i = 0; i<dim_sim;++i)
    {
        sim[i] = new double[dim];
    }
}

simplex::~simplex()
{
    for( int i = 0; i < dim_sim; ++i)
    {
        delete[] sim[i];
    }
    delete[] sim;
}

void simplex::update(double** points)
{
    sim=points;
}

double** simplex::getSimp()
{
    return sim;
}

int simplex::getDim()
{
    return dim;
}

int simplex::getDimS()
{
    return dim_sim;
}

void simplex::init(double* starting_point)
{   
    //inicjalizacja simpleksu
    for(int i = 0; i < dim; ++i)
    {
        sim[0][i]=starting_point[i];
    }

    //uzupelnienie o inne punkty
    for(int j = 1; j < dim_sim; ++j)
    {
        for(int i = 0, k = 1; i < dim; ++i, ++k)
        {
            sim[j][i]=sim[0][i];
            sim[k][i]=sim[0][i]+0.1;
        }
    }
    return;
}

double calculate_mean(const double* f_values, int dim)
{
    double mean = 0;
    for(int i = 0; i < dim; ++i)
    {
        mean += f_values[i];
    }
    return mean/dim;
}

double* calculate_mean_2d(double** simp, int f_dim)
{
    double* mean =  new double[f_dim];

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

double calculate_std(const double* f_values, int dim)
{
    double mean = calculate_mean(f_values, dim);
    double variance = 0;
    for( int i = 0; i < dim; ++i)
    {
        variance += pow((f_values[i] - mean), 2);
    }
    return sqrt(variance/(dim-1));
}

void idx_sort(double* f_values, simplex& simp)
{
    int dim_s = simp.getDimS();    // liczba punktów w simpleksie
    int dim = simp.getDim();       // wymiar przestrzeni
    double** sim = simp.getSimp(); // wskaźnik do macierzy punktów (simpleksu)

    // Utwórz wektor indeksów
    std::vector<int> idx(dim_s);
    for (int i = 0; i < dim_s; ++i)
    {
        idx[i] = i;
    }
    // Posortuj indeksy względem f_values
    std::sort(idx.begin(), idx.end(), [&](int i, int j) {
        return f_values[i] < f_values[j];
    });

    // Stwórz kopie posortowanych danych
    double* f_sorted = new double[dim_s];
    double** sim_sorted = new double*[dim_s];
    for (int i = 0; i < dim_s; ++i) {
        f_sorted[i] = f_values[idx[i]];

        sim_sorted[i] = new double[dim];
        for (int j = 0; j < dim; ++j) {
            sim_sorted[i][j] = sim[idx[i]][j];
        }
    }

    // Zastąp dane oryginalne
    for (int i = 0; i < dim_s; ++i) {
        f_values[i] = f_sorted[i];
        for (int j = 0; j < dim; ++j) {
            sim[i][j] = sim_sorted[i][j];
        }
        delete[] sim_sorted[i];
    }
    delete[] sim_sorted;
    delete[] f_sorted;
}

void printSimp(simplex* simp, funkcja func)
{
    int simplex_point_no = simp->getDimS();
    int dimension_no = simp->getDim();

    double** pts = simp->getSimp();
    double f;

    for(int i = 0; i < simplex_point_no; ++i)
    {
        for( int j = 0; j < dimension_no; ++j)
        {
            printf("%6.3f ", pts[i][j]);
        }
        f=func(pts[i]);
        printf("%6.3f \n", f);
    }
}

double* min_search(const int no_dim, double* start_pt, funkcja function)
{
    simplex starting_point(no_dim);
    int simplex_point_no = starting_point.getDimS();
    int dimension_no = starting_point.getDim();
   
    starting_point.init(start_pt);
    //printSimp(&starting_point, function);
    double** current_simplex = starting_point.getSimp();

    //inicjalizacja wektora funkcji
    double* f_values = new double[starting_point.getDimS()];
    for(int i = 0; i < simplex_point_no; ++i)
    {
        f_values[i]=function(current_simplex[i]);
    }

    //zmienne pomocnicze
    double* x_o = new double [dimension_no];
    double* x_r = new double [dimension_no];
    double* x_c = new double [dimension_no];
    double* x_e = new double [dimension_no];
    double f_r, f_c, f_e;

    //glowna pętla
    for(int i = 0; i < MAX_I; ++i)
    {
        //warunek wyjscia (nie dziala)
        if(calculate_std(f_values, simplex_point_no) < EPSILON) 
        {
            //printf("%s", "osiagnieto pozadane odchylenie standradowe \n");
            break;
        }

        idx_sort(f_values, starting_point);

        x_o = calculate_mean_2d(current_simplex, dimension_no);

        for(int j = 0; j < dimension_no; ++j)
        {
            x_r[j] = x_o[j] + ALPHA * (x_o[j] - current_simplex[simplex_point_no-1][j]);
        }
        f_r=function(x_r);

        if(f_r >= f_values[0] && f_r < f_values[simplex_point_no-1])
        {
            for (int j = 0; j < dimension_no; ++j)
                current_simplex[simplex_point_no-1][j] = x_r[j];
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
                    current_simplex[simplex_point_no-1][j] = x_e[j];
                f_values[simplex_point_no-1] = f_e;
                continue;
            }
            else
            {
                for (int j = 0; j < dimension_no; ++j)
                    current_simplex[simplex_point_no-1][j] = x_r[j];
                f_values[simplex_point_no-1] = f_r;
                continue;
            }
        }
        else if(f_r >= f_values[simplex_point_no-1])
        {
            for(int j = 0; j < dimension_no; ++j)
            {
                x_c[j] = x_o[j] + BETA * (current_simplex[simplex_point_no-1][j] - x_o[j]);
            }
            f_c=function(x_c);
            
            if(f_c < f_r)
            {
                for (int j = 0; j < dimension_no; ++j)
                    current_simplex[simplex_point_no-1][j] = x_c[j];
                f_values[simplex_point_no-1] = f_c;
                continue;
            }
        }

        for(int j = 1; j < simplex_point_no; ++j)
        {   
            for(int k = 0; k < dimension_no; ++k)
            {
                current_simplex[j][k] = current_simplex[0][k] + DELTA * (current_simplex[j][k] - current_simplex[0][k]);
            }
            f_values[j] = function(current_simplex[j]);
        }
    }

    delete[] x_o;
    delete[] x_r;
    delete[] x_c;
    delete[] x_e;
    delete[] f_values;
    //printSimp(&starting_point, function);
    return current_simplex[0];
}