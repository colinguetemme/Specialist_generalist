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
int L; //period
double A; //amplitude of fluctuation
double d; //magnitude of stochasticity

//////////////////
//// VARIABLES ///
//////////////////

vector<double> opti_G; // Simulation time  
vector<double> r; // Recombination rates
const int nn = n - 1; // Number of loci -1
int ng = pow(2, n);// Number of gametes
int ngg = pow(2, nn);// Number of recombination positions
double opt; // Optimum
double alpha[n]; // Allelic values
vector<double> p; // Gamete frequencies

struct init{
	vector<double> loci_values;	// The positive alleles values of the population
    vector<double> gamete_values; // Phenotypic values of gametes
	vector<double> gamete_distr; // Initial distribution of gametes
	vector<string> gamete_scheme; // Binary of each gametes
	vector<vector<vector<double>>> rec_table; // 3D table of every recombination probability
};

struct one_step{
	vector<double> new_distr; // The new distribution from eqn 4 Burger
	double mean_fitness; // The mean fitness of the population at that step
	double var_genetics; // The genetic variance of the population at that step
};


// A set of step of the simulation (more likely the full set of the simulation)
struct many_steps{
	vector<vector<double>> all_distr; // The new distribution from eqn 4 Burger
	vector<double> mean_fitness; // The mean fitness of the population at each time step
	vector<double> gamete_values; // The gamete values of the population
	vector<double> loci_values; // The positive allele values of the population 
};

// Structure with the statistical results as in the article of Burger(2002)
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
					   vector<double>, double, vector<string>, double, double, bool);
many_steps one_simul(double, double, int, int, double, vector<double>, bool);
init initialisation(bool);
vector<string> permut_gamete(string j, string k, vector<int>);
double jk_i_recombination(string, string, string, vector<double>);
vector<double> environment(double, double, int, bool);
vector<double> environment2(vector<double> , int, double, bool);
results statistics(vector<vector<double>>, vector<double>, vector <double>, vector<double>, int);
