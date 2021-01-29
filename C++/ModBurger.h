#pragma once

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
#include <bits/stdc++.h> 

using namespace std; 

//Parameters
int t; //simulation time  

const int n = 4;
const int nn = n - 1;
int ng = pow(2, n);//number of gametes

double u; //mutation rate
double s; //strength of selection
vector<double> r; //recombination rates

int L; //period
double A; //amplitude of fluctuation
double d; //magnitude of stochasticity

//Variables
double opt; //optimum


double alpha[n]; //allelic values

vector<double> p; //gamete frequencies

//table of recombination probabilities
double*** rec_table;

//struct rec_table {
//	double jk_table[n][n];  
//};
//vector<rec_table> rec_tables; //gamete values


vector<std::string> gametes;
vector<double> genotypes;


// Random numbers generators
// seed random number generator
std::random_device rd;
std::mt19937 rdgen(rd());

std::uniform_real_distribution<> unif(0.0, 1.0);


//Functions
void initialisation(void);
double recombination(string, string, string);
vector<string> permut_gamete(string j, string k, vector<int>);
double jk_i_recombination(string, string, string, vector<double>);
const string Int2Str(const int x);
void rec_test(void);

