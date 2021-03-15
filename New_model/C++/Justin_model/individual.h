
#pragma once

# include "chromosome.h"

using namespace std;

// Maybe will remove since we do not need diploid for now //!//
class individual{
    public:
chromosome chr;

double init_fit_val = 1; // the initial fitness value
double fit_val; // the current fitness value

void fit_value(double);
};

// The class containing a vector of individual
class population{
    public:
    population(int);
    ~population();
vector<individual> ind;
void new_generation(double, para_ind);
};