//////////////
// INCLUDES //
//////////////

#include "ModBurger.h"

///////////////////
// MAIN FUNCTION //
///////////////////

int main(){

	////////////////
	// PARAMETERS //
	////////////////

	// a is the amplitude of the sinusoid
	// L is the period of the sinusoid
	// d is the stochastic parameter
	// s is the strength of the selection
	// u is the mutation rate
	// n is the number of loci 

	results result;
	double a = 0.5;
	int L = 50;
	double d = 0;
	double s = 5;
	double u = 0 ; //pow(10, -6);
	const int n = 2;

	////////////////////
	
	int nn = n - 1;
	int ng = pow(2, n); //number of gametes
	int ngg = pow(2, nn); // number of recombination positions

	/////////////////
	// OUTPUT FILE //
	/////////////////

	// Setting up the output file
	// Use the R code to quick plot the results

	string filename = "../results/full_distr/";
	
	filename.append("a");
	filename.append(to_string(int(a*100)));
	filename.append("_L");
	filename.append(to_string(L));
	filename.append("_d");
	filename.append(to_string(int(d*100)));
	filename.append("_s");
	filename.append(to_string(int(s)));
	filename.append("_u");
	if (u == 0){
		filename.append(to_string(0));
	} else {
		filename.append(to_string(int(log10(u))));
	}
	filename.append("_n");
	filename.append(to_string(n));
	filename.append(".csv");	
	std::ofstream outfile (filename);

	outfile << "sim_num" << "\t";
	outfile << "time" << "\t";
	outfile  << "gamete_num" << "\t";
	outfile  << "prop" << "\t";
	outfile  << "gamete_value" << "\n";


	string filename_stat = "../results/stats/";

	filename_stat.append("a");
	filename_stat.append(to_string(int(a*100)));
	filename_stat.append("_L");
	filename_stat.append(to_string(L));
	filename_stat.append("_d");
	filename_stat.append(to_string(int(d*100)));
	filename_stat.append("_s");
	filename_stat.append(to_string(int(s)));
	filename_stat.append("_u");
	if (u == 0){
		filename_stat.append(to_string(0));
	} else {
		filename_stat.append(to_string(int(log10(u))));
	}
	filename_stat.append("_n");
	filename_stat.append(to_string(n));
	filename_stat.append(".csv");
	std::ofstream outfile_stat (filename_stat);

	outfile_stat << "sim_num" << "\t";
	outfile_stat << "stat" << "\t";
	outfile_stat << "value" << "\n";

	/////////////////
	// SIMULATIONS //
	/////////////////

	// The number of simulation for one set of parameter
	// is suppose to be 4000 in the paper of Burger (//!// very long simulation time)

	int num_sim = 100;
	vector<vector<double>> (4, vector<double> (num_sim, 0));

	for (int e = 0; e<num_sim; e++){
		data = many_steps();
		cout << "\r Processing Simulations: " << e+1 << "/" << num_sim;

		// Creates the varying environment
		vector<double> opti_G = environment(a, d, L);

		// Do the simulation for this environment
		data = one_simul(s, u, n, L, d, opti_G);
		
		int keep_end = 0;
		if (d>0){
			keep_end = 50000 - (10*L);
		}

		for (int i = 0; i<data.mean_fitness.size(); i++){
			outfile << e << "\t" << i << "\t" << -1 << "\t" << opti_G[i+keep_end] <<  "\t" << 0.1 << "\n";
		}

		for (int i = 0; i<data.mean_fitness.size(); i++){
			for (int j = 0; j<ng; j++){
				outfile << e << "\t" << i << "\t" << j << "\t" << data.all_distr[i][j] << "\t" << data.gamete_values[j] << "\n";
			}
		}

		// Calculates the values of interest
		result = statistics(data.all_distr, data.gamete_values, data.loci_values, data.mean_fitness, n);

		// Add the results to the ouput file
		outfile_stat << e << "\t" << "mean_gen" << "\t" << result.mean_gen << "\n";
		outfile_stat << e << "\t" << "var_gen" << "\t" << result.var_gen << "\n";
		outfile_stat << e << "\t" << "ratio" << "\t" << result.ratio << "\n";
		outfile_stat << e << "\t" << "geom_fitness" << "\t" << result.geom_fitness << "\n";
		
	}

	outfile.close();
	outfile_stat.close();

	cout << "\n\n" << "DONE" << "\n" << endl;
	cout << "Find the results at: " << filename_stat << "\n" << "\n" << endl;
	
	
	return 0;
}

///////////////////
// ENVIRONMENT() //
///////////////////
// Inputs:
//	'a' the amplitude of the sinusoid
//	'd'	the stochastic factor
//	'L' the period (number of generations)
//
// Outputs:
//	'opti_G' a vector of double containing the value
//		optimum genotype through time 


vector<double> environment(double a, double d, int L)
{	
	int tmax;
	if (d==0){
		// when d == 0, the environment is cyclic with period L 
		// so we just need the values in one cycle
		tmax = L;
	} else {
		// when d > 0, the simulation stop at 50000 generations Burger(2002)
		tmax = 50000;
	}

	vector<double> opti_G(tmax, 0);
	double var = d * a;

	// a rng following a normal distribution
	std::normal_distribution<> norm(0, var);

	for (int t = 0; t < tmax; t++)
	{
		opti_G[t] = norm(rdgen) + 0.5 + a * sin(2 * M_PI * t / L);
	}

	return opti_G;
}

/////////////////
// ONE_SIMUL() //
/////////////////
// Inputs:
//	's' the strength of selection
//	'u'	the mutation rate
//	'n' the number of loci
//  'opti_G' the optimum genotype
//
// Outputs:
//	'data'  //!//

many_steps one_simul(double s, double u, int n, int L, double d, vector<double> opti_G){
	
	init init_values = initialisation();

	// Initialze the variable which will stock the distributions
	// of the gametes through the time
	one_step results;
	many_steps data;

	
	// If there is no stochasticity (meaning d==0)
	if (d==0){

		vector<vector<double>> all_distr(L, vector<double>(ng, 0));
		vector<double> mean_fit(L, 0);


		// initialize the distribution
		data.all_distr = all_distr;
		data.mean_fitness = mean_fit;
		data.gamete_values = init_values.gamete_values;
		data.loci_values = init_values.loci_values;

		vector<double> old_distr = init_values.gamete_distr;


		// Does the simulation
		for (int t = 0; t<300000; t++){
			results = new_distributions(init_values.gamete_values, init_values.rec_table, old_distr, u, init_values.gamete_scheme, opti_G[t%L], s);

			old_distr = results.new_distr;

			// test whenever we need to check the stability
			if ((t%L == 0) & (L>0)){
	
				// Calculate the distance;
				double distance = 0;
				for (int i = 0; i<ng; i++){
					distance += pow(old_distr[i] - data.all_distr[0][i], 2);
				}
				distance = sqrt(distance);

				// Check if the distance is higher than 10^-12
				if (distance < pow(10, -12)){

					return data;
				}
			}

			// Update the distributions
			data.all_distr[t%L] = results.new_distr;
			data.mean_fitness[t%L] = results.mean_fitness;

			double sum = 0;
			for (int i = 0; i<ng; i++){
				sum += data.all_distr[t%L][i];
			}
	
		}
	}

	// if there is stochasticity
	if (d>0){
		vector<vector<double>> all_distr(10*L, vector<double>(ng, 0));
		vector<double> mean_fit(10*L, 0);


		// initialize the distribution
		data.all_distr = all_distr;
		data.mean_fitness = mean_fit;
		data.gamete_values = init_values.gamete_values;
		data.loci_values = init_values.loci_values;


		vector<double> old_distr = init_values.gamete_distr;

		// Update the distributions
		for (int t = 0; t<50000; t++){
			results = new_distributions(init_values.gamete_values, init_values.rec_table, old_distr, u, init_values.gamete_scheme, opti_G[t], s);
			data.all_distr[t%(L*10)] = results.new_distr;
			data.mean_fitness[t%(L*10)] = results.mean_fitness;
			old_distr = results.new_distr; 
		}
	}

	return data;
}

//////////////////
// STATISTICS() //
//////////////////


// Inputs:
//	'gamete_distr' the 2D vector containing the distribution 
// 		of each gametes for the kept generations
//	'gamete_values'	the genotypic value of each gamete
//	'loci_value' the genotypic value of the positive allele
//		at each loci
//	'mean_fitness' the mean fitness of the population for the
//		kept generations
//
// Outputs:
//	'result.':
//		- 'mean_gen' the averaged value of the mean genotype in the population
//		- 'var_gen' the averaged genetic variance
//		- 'ratio' the ratio between the variance and Vmax (the maximum possible variance)
//		- 'geom_mean' the geometric mean of the fitness


results statistics(vector<vector<double>> gamete_distr,
 vector<double> gamete_values,
 vector <double> loci_values,
  vector<double> mean_fitness,
  int n){

	results result;

	// The averaged value of the mean genotype in the population
	double mean_gen = 0;
	vector<double> mean_gen_gen(gamete_distr.size(), 0);
	for (int i = 0; i<gamete_distr.size(); i++){
	    for (int j = 0; j<ng; j++){
		    mean_gen += gamete_distr[i][j] * gamete_values[j];
			mean_gen_gen[i] += gamete_distr[i][j] * gamete_values[j];
	    }
	}
	
	mean_gen /= gamete_distr.size();

	result.mean_gen = mean_gen;

	// The averaged genetic variance
	// Need to check if var for each generation or for each cycle //!//
	double var_gen = 0;
	for (int i = 0; i<gamete_distr.size(); i++){
	    for (int j = 0; j<ng; j++){

	        var_gen += pow(gamete_distr[i][j] * gamete_values[j] - mean_gen_gen[i], 2);
	    }
	}
	var_gen /= gamete_distr.size();
	result.var_gen = var_gen;

	// The ratio between the variance and Vmax (the maximum possible variance)
	double Vmax = 0;
	for (int i = 0; i<n; i++){
		Vmax += pow(loci_values[i], 2);
	}
    double ratio = var_gen / (Vmax/2) / ng;
	result.ratio = ratio;

	// The geometric mean of the fitness
	double geom_fitness = 1;
	double root = gamete_distr.size();
	for (int i = 0; i<gamete_distr.size(); i++){
		geom_fitness *= pow(mean_fitness[i], 1/root);
	}
	result.geom_fitness = geom_fitness;

	return result;
}

//////////////////////
// INITIALISATION() //
//////////////////////


// Inputs:
//
// Outputs:
//	'init_values.':
//		- 'loci_values' the genotypic value of the positive allele at each loci
//		- 'gamete_values' the genotypic value of each gamete
//		- 'gamete_scheme' the allelic composition of each possible gamete (with
//				0, the allele of value 0, and 1, the positive allele)
//		- 'rec_table' the 3D vector containting the probability to obtain gamete i
//				from gamete j and k by recmobination


init initialisation(void)
{	
	init init_values;
	double sum = 0.0;


	vector<double> alpha(n);
	alpha[0] = unif(rdgen);

	sum += alpha[0];

	for (int i = 0; i < n - 1; i++)
	{
		alpha[i + 1] = unif(rdgen); //initialise alleles
		sum += alpha[i + 1];
		r.push_back(0.0);
		r[i] = unif(rdgen) / 2.0; //initialise recombination rate
	}

	for (int i = 0; i < n; i++)
	{
		alpha[i] = 0.5 * alpha[i] / sum; //scale alleles

	}

	vector<string> gametes;
	vector<double> genotypes;

	//initialise gametes
	for (int i = 0; i < ng; i++)
	{

		std::bitset<n> alleles(i);
		gametes.push_back(alleles.to_string<char, std::string::traits_type, std::string::allocator_type>());
		sum = 0.0;

		for (int j = 0; j < n; j++)
		{
			if (((int)gametes[i][j]) - 48 == 1)
				sum += alpha[j];
		}
		genotypes.push_back(sum);
	}

	//initiliase recombination table
	vector<vector<vector<double>>> rec_table(ng, (vector<vector<double>>(ng, (vector<double>(ng, 0)))));
	for (int j = 0; j < ng; j++)
	{
		for (int k = 0; k < ng; k++)
		{
			for (int i = 0; i < ng; i++)
			{
				rec_table[j][k][i] = jk_i_recombination(gametes[i], gametes[j], gametes[k], r);
			}
		}
	}

	sum = 0;
	vector<double> distr;
	double x;
	for (int i = 0; i<ng; i++){
		x = unif(rdgen);
		distr.push_back(x);
		sum += x;
	}
	for (int i = 0; i<ng; i++){
		distr[i] /= sum;
	}

	init_values.loci_values = alpha;
	init_values.gamete_values = genotypes;
	init_values.gamete_distr = distr;
	init_values.gamete_scheme = gametes;
	init_values.rec_table = rec_table;
	
	return init_values;
}

//////////////////////////
// JK_I_RECOMBINATION() //
//////////////////////////
// Inputs:
//	'i' the string corresponding to gamete i
//	'j' the string corresponding to gamete j
//	'k' the string corresponding to gamete k
//	'r' the vector of rate of recombination
//
// Outputs:
//	'R' the probability to goes from jk to i by recombination

double jk_i_recombination(string i, string j, string k, vector<double> r)
{

	double R = 0.0; // final probability of recombination btw j, k, and i
	double prod = 1.0;
	//int ngg = pow(2.0, nn); // number of possible recombinations combinations
	int **t; // table of all recombinations

	vector<int> rec_vec; // all recombinations possibilities
	vector<string> jk; // binary of new j and k

	t = new int *[ngg];
	for (int z = 0; z < ngg; z++)
	{ // z go trough all cmb of recombinations
		rec_vec.clear();
		t[z] = new int[nn];
		std::bitset<nn> allele(z);
		for (int zz = 0; zz < nn; zz++)
		{ // zz go trhough all the rec_pos
			t[z][zz] = allele[zz] * (zz + 1);
			if (t[z][zz] > 0)
				rec_vec.push_back(t[z][zz]);
		}

		if (z > 0)
		{
			jk = permut_gamete(j, k, rec_vec);
		}
		else
		{
			jk.push_back(j);
			jk.push_back(k);
			for (int zz = 1; zz < nn; zz++)
				rec_vec.push_back(-zz);
		}

		for (int a = 0; a < 2; a++)
		{
			if (jk[a] == i)
			{

				int rv = 0;
				prod = 1.0;
				for (int zz = 0; zz < (int)r.size(); zz++)
				{
					if (zz == rec_vec[rv] - 1)
					{
						prod *= r[zz];
						rv++;
					}
					else
						prod *= 1.0 - r[zz];
				}
				R += 0.5 * prod;
			}
		}
	}
	return R;
}

//////////////////////////
// PERMUT_GAMETE() //
//////////////////////////
// Inputs:
//	'i' the string corresponding to gamete i
//	'j' the string corresponding to gamete j
//	'k' the string corresponding to gamete k
//	'r' the vector of rate of recombination
//
// Outputs:
//	'R' the probability to goes from jk to i by recombination

vector<string> permut_gamete(string j, string k, vector<int> rec)
{
	int rec_spot;

	string j_bin_rec, k_bin_rec;
	vector<string> jk;

	rec_spot = rec[0];

	for (int z = 0; z < rec_spot; z++)
		j_bin_rec += j[z];
	for (int z = rec_spot; z < (int)j.size(); z++)
		j_bin_rec += k[z];

	for (int z = 0; z < rec_spot; z++)
		k_bin_rec += k[z];
	for (int z = rec_spot; z < (int)k.size(); z++)
		k_bin_rec += j[z];

	rec.erase(rec.begin());

	if ((int)rec.size() > 1)
		permut_gamete(j_bin_rec, k_bin_rec, rec);

	jk.push_back(j_bin_rec);
	jk.push_back(k_bin_rec);

	return jk;
}

one_step new_distributions(vector<double> gamete_values, vector<vector<vector<double>>> dtf_rec,
					   vector<double> gamete_distr, double mut_rate, vector<string> gam_bin, double opti_Gt, double strength)
{
	int ng = pow(2, n);
	vector<double> p_star(ng, 0);

	double sum = 0;
		for (int i = 0; i<ng; i++){
			sum += gamete_distr[i];
		}

	double mean_fitness = 0;
	for (int j = 0; j < ng; j++)
	{
		for (int k = j; k < ng; k++)
		{

			mean_fitness += exp(-strength * pow(gamete_values[j] + gamete_values[k] - opti_Gt, 2)) * gamete_distr[j] * gamete_distr[k];
		}
	}
	
	for (int i = 0; i < ng; i++)
	{
		for (int j = 0; j < ng; j++)
		{
			for (int k = j; k < ng; k++)
			{
				p_star[i] += exp(-strength * pow(gamete_values[j] + gamete_values[k] - opti_Gt, 2)) * gamete_distr[j] * gamete_distr[k] * dtf_rec[i][j][k] / mean_fitness;
			}
		}
	}

	double mutation_ij = 0;
	
	for (int i = 0; i < ng; i++)
	{
		double tot_mut = 0;
		for (int j = 0; j < ng; j++)
		{
			
			int nm = 0;
			if (i != j)
			{
				for (int z = 0; z < n; z++)
				{
					if (gam_bin[i][z] != gam_bin[j][z]) nm++;
				}
				mutation_ij = pow(mut_rate, nm);
				tot_mut += tot_mut + p_star[j] * mutation_ij - p_star[i] * mutation_ij;
			}
		}
		p_star[i] += tot_mut;
	}

	// TEST NOT SURE
	double sum2 = 0;
	for (int i = 0; i<ng; i++){
		sum2+= p_star[i];
	}
	
	for (int i = 0; i<ng; i++){
		p_star[i] /= sum2;
	}
	
	double var = 0;

	
	one_step next_step;
	next_step.new_distr = p_star;
	next_step.mean_fitness = mean_fitness;
	
	next_step.var_genetics = var;
	return next_step;
}

