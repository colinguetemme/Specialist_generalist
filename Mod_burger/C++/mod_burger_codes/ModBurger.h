#pragma once

/////////////////
/// LIBRARIES ///
/////////////////

#include <stdio.h>
#include <stdlib.h>
//#include <io.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <numeric>
#include <time.h>
#include <random>
#include <iterator>
#include <bitset> 

using namespace std; 

//////////////////
/// PARAMETERS ///
//////////////////

/// Genetics
const int n = 4; // number of loci
double u; //mutation rate
double s; //strength of selection

/// Temporal variation
int tmax;
int L; //period
double A; //amplitude of fluctuation
double d; //magnitude of stochasticity

//////////////////
//// VARIABLES ///
//////////////////

vector<double> opti_G; //simulation time  
vector<double> r; //recombination rates
const int nn = n - 1;
int ng = pow(2, n);//number of gametes
int ngg = pow(2, nn);// number of recombination positions
double opt; //optimum
double alpha[n]; //allelic values
vector<double> p; //gamete frequencies

struct init{
	vector<double> loci_values;
    vector<double> gamete_values; // phenotypic values of gametes
	vector<double> gamete_distr; // initial distribution of gametes
	vector<string> gamete_scheme; // binary of each gametes
	vector<vector<vector<double>>> rec_table; // 3d table of every recombination probability
};

struct one_step{
	vector<double> new_distr; // the new distribution from eqn 4 Burger
	double mean_fitness;
	double var_genetics;
};

struct many_steps{
	vector<vector<double>> all_distr; // the new distribution from eqn 4 Burger
	vector<double> mean_fitness;
	vector<double> gamete_values;
	vector<double> loci_values;
};

struct results{
	double mean_gen;
	double var_gen;
	double ratio;
	double geom_fitness;
};

many_steps data;
vector<std::string> gametes;
vector<double> genotypes;

// Random numbers generators
// seed random number generator
std::random_device rd;
std::mt19937 rdgen(rd());
std::uniform_real_distribution<> unif(0.0, 1.0);


/////////////////
/// FUNCTIONS ///
/////////////////

one_step new_distributions(vector<double>, vector<vector<vector<double>>>,
					   vector<double>, double, vector<string>, double, double);
many_steps one_simul(double, double, int, int, double, vector<double>);
init initialisation();
vector<string> permut_gamete(string j, string k, vector<int>);
double jk_i_recombination(string, string, string, vector<double>);
vector<double> environment(double, double, int);
results statistics(vector<vector<double>>, vector<double>, vector <double>, vector<double>, int);
