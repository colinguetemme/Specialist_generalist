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
double u = 0.00005; //mutation rate
double s = 10; //strength of selection

/// Temporal variation
int tmax = 200;
int L = 24; //period
double A = 0.5; //amplitude of fluctuation
double d = 0.1; //magnitude of stochasticity

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
    vector<double> gamete_values; // phenotypic values of gametes
	vector<double> gamete_distr; // initial distribution of gametes
	vector<string> gamete_scheme; // binary of each gametes
	vector<vector<vector<double>>> rec_table; // 3d table of every recombination probability
};


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

vector<double> new_distributions(vector<double>, vector<vector<vector<double>>>,
					   vector<double>, double, vector<string>, double, double);
init initialisation();
vector<string> permut_gamete(string j, string k, vector<int>);
double jk_i_recombination(string, string, string, vector<double>);
vector<double> environment(int, int, double, double);
