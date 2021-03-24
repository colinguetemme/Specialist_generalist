/**
 * @file ModBurger.cpp
 *
 * @author Colin Guetemme (colin.guetemme@live.fr)
 * 
 * @brief Do the simulation for the same model as the one from Burger(2002), for the purpose of 
 * our research, also available for haploid population (i.e. in this model, no recombination).
 * Briefly, the model test the genetic variation of a population evolving in a time-varying environment.
 * The hypothesis being that a high rate of variation will leads to a generalist species and therefore
 * very high genetic variance and low rate to a specialist species, because this adapt has the time to 
 * specialise to the new environment each time.
 * 
 * @version 0.1
 * @date 2021-03-22
 * @copyright Copyright (c) 2021
 * 
 */

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

	results result;
	double a = 0.5; // Amplitude of the sinusoid
	int L = 50;	// Period of the sinusoid
	double d = 0; // Stochasticity coefficient
	double s = 5; // Strength of selection
	double u = 0 ; // Mutation rate
	// double u = pow(10, -6); 
	const int n = 2; // Number of loci
	bool diploid = 0; // If the population is haploid (= 0) or diploid (= 1)

	////////////////////
	
	int nn = n - 1;
	int ng = pow(2, n); // Number of gametes
	int ngg = pow(2, nn); // Number of recombination positions

	//////////////////
	// OUTPUT FILES //
	//////////////////

	// Setting up the output file
	// The output have been organize as a ggplot2 data frame (in R)
	// Use the R code to quick plot the results

	string file;

	// Creates the file name after the parameters of the simulation

	if (diploid == 1){
		file = "dipl";
	} else {
		file = "hapl";
	}
	file += "_a" + to_string(int(a*100)) + "_L" + to_string(L) + "_d" + to_string(int(d*100)) +
	 "_s" + to_string(int(s)) + "_n" + to_string(n) + "_u";

	if (u == 0){
		file += to_string(0);
	} else {
		file += to_string(int(log10(u)));
	}
	
	file += ".csv";	

	// CREATE THE FILES

	// The file stored in full_distr will contain the distribution of each gamete
	// at every time step ketp for the results
	string dirname = "../results/full_distr/";
	string filename = dirname + file;

	// The file stored in stats will contain the statistical results of the simulation
	string stat_dirname = "../results/stats/";
	string filename_stat = stat_dirname + file;

	std::ofstream outfile (filename);
	std::ofstream outfile_stat (filename_stat);
	
	// Headers of the csv
	outfile_stat << "sim_num" << ",";
	outfile_stat << "stat" << ",";
	outfile_stat << "value" << "\n";

	/////////////////
	// SIMULATIONS //
	/////////////////

	// The number of simulation for one set of parameter
	// is suppose to be 4000 in the paper of Burger (//!// very long simulation time)
	int num_sim = 100;

	for (int e = 0; e<num_sim; e++){
		data = many_steps(); // gather the results of one simulaion
		cout << "\r Processing Simulations: " << e+1 << "/" << num_sim;

		// The optimum genotype in a varying environment 
		vector<double> opti_G = environment(a, d, L, diploid);

		// Do the simulation for this environment
		data = one_simul(s, u, n, L, d, opti_G, diploid);
		
		int keep_end = 0; // Only important for stochastic, index for the last values of the environment
		if (d>0){
			keep_end = 50000 - (10*L); // Start from the last 10 cycles
		}

		// In the csv, the environment will be added to same dataframe as the other gamete distribution
		// but either of having a gamete index it will have the index -1
		for (int i = 0; i<data.mean_fitness.size(); i++){
			outfile << e << "," << i << "," << -1 << "," << opti_G[i+keep_end] <<  "," << 0.1 << "\n";
		}

		// Add the gamete distribution for each gamete
		for (int i = 0; i<data.mean_fitness.size(); i++){
			for (int j = 0; j<ng; j++){
				outfile << e << "," << i << "," << j << "," << data.all_distr[i][j] << "," << data.gamete_values[j] << "\n";
			}
		}

		// Calculates the values of interest
		result = statistics(data.all_distr, data.gamete_values, data.loci_values, data.mean_fitness, n);

		// Add the statistics results to the ouput file
		outfile_stat << e << "," << "mean_gen" << "," << result.mean_gen << "\n";
		outfile_stat << e << "," << "var_gen" << "," << result.var_gen << "\n";
		outfile_stat << e << "," << "ratio" << "," << result.ratio << "\n";
		outfile_stat << e << "," << "geom_fitness" << "," << result.geom_fitness << "\n";
		
	}

	// // // // END OF THE SIMULATION // // // // 

	outfile.close();
	outfile_stat.close();

	cout << "\n\n" << "DONE" << "\n" << endl;
	cout << "Find the results at: " << filename_stat << "\n" << "\n" << endl;
	
	
	return 0;
}


// // // OTHER FUNCTIONS // // //

/////////////////
// ENVIRONMENT //
/////////////////

/**
 * @brief Creates a fluctuating environment based on a sinusoid function,
 * 
 * if there is no stochasticity (d=0), the output will corresponds to one cycle,
 * if there is stochasticity (d>0), the output will corresponds to 50000 generations.
 * Also for the stocastic environment, the optimum is sampled from a normal distribution
 * N(opti, d*a), where opti is the optimum without stochasticity.
 * 
 * @param a the amplitude of the sinusoid
 * @param d the stochasticity coefficient
 * @param L the period (number of generation) of the sinusoid
 * @return a vector with the optimum genotype at each time 
 */
vector<double> environment(double a, double d, int L, bool dipl)
{	
	int tmax; // The maximum number of generation to look at
	if (d==0){
		// when d == 0, the environment is cyclic with period L 
		// so we just need the values in one cycle
		tmax = L;
	} else {
		// when d > 0, the simulation stop at 50000 generations, Burger(2002)
		tmax = 50000;
	}

	vector<double> opti_G(tmax, 0);
	double var = d * a;

	// a rng following a normal distribution
	std::normal_distribution<> norm(0, var);

	double coeff; // if diploid one gamete should count for 0.5 of the genotypic value;
	if (dipl == 1){
		coeff = 0.5;
	} else {
		coeff = 1;
	}
	for (int t = 0; t < tmax; t++)
	{
		opti_G[t] = norm(rdgen) + coeff + a * sin(2 * M_PI * t / L);
	}

	return opti_G;
}

/**
 * @brief More flexible environment, not perfect but gives the possibility to have non sinus function.
 * //!// The period should be changed to have a period rather than a time in each state.
 * For the stochasticity get 50000 steps and for the non-stochastic just generate one cycle (L*states.size()).
 * 
 * @param states vector containing the states of the changing environment
 * @param L in this case, corresponds to the time spend in each state, the real period being L * n_states
 * @param d the stochastic parameters, 
 * @param dipl population is diploid (= 1) or haploid (= 0)
 * @return The optimum genotype at each time step
 */
vector<double> environment2(vector<double> states, int L, double d, bool dipl){ 
	if (d>0){
		
		vector<double> opti_G(0, 50000);

		double mean_states = 0;
		double var_states = 0;
		for(int i = 0; i>states.size(); i++){
			mean_states += states[i];
		}
		mean_states /= states.size();
		for(int i = 0; i>states.size(); i++){
			var_states += pow((states[i]-mean_states),2);
		}

		std::normal_distribution<> norm(0, var_states);

		int index;
		for (int i = 0; i<50000; i++){
			index = floor(i%(states.size()*L)/L); // The time step in a cycle
			opti_G[i] = states[index] + norm(rdgen);
		}
	
		return opti_G;

	} else {
		int index;
		for (int i = 0; i>(L*states.size()); i++){
			index = floor(i%(states.size()*L)/L); // The time step in a cycle
			opti_G[i] = states[index];
		}
	}
}

///////////////
// ONE_SIMUL //
///////////////

/**
 * @brief does one simulation of the model with the given parameters
 * 
 * @param s the strength of selection
 * @param u the mutation rate
 * @param n the number of loci
 * @param L the number of generation for one cycle (period)
 * @param d the stochasticity coefficient of the environment
 * @param opti_G the vector of optimum genotype through time
 * 
 * @return a structure containing:
 * 		- a 2D vector of the distribution of each gamete at each time step
 * 		- a vector of the mean fitness at each timestep
 * 		- a vector of the gamete values
 * 		- a vector of the loci values
 */
many_steps one_simul(double s, double u, int n, int L, double d, vector<double> opti_G, bool dipl){
	
	init init_values = initialisation(dipl);

	// Initialze the variable which will stock the distributions
	// of the gametes through the time
	one_step results;
	many_steps data;

	
	// WITHOUT STOCHASTICITY (meaning d==0)

	if (d==0){

		vector<vector<double>> all_distr(L, vector<double>(ng, 0));
		vector<double> mean_fit(L, 0);


		// initialize the distribution
		data.all_distr = all_distr;
		data.mean_fitness = mean_fit;
		data.gamete_values = init_values.gamete_values;
		data.loci_values = init_values.loci_values;

		vector<double> old_distr = init_values.gamete_distr;


		// DO THE SIMULATION //

		for (int t = 0; t<300000; t++){
			results = new_distributions(init_values.gamete_values, init_values.rec_table, old_distr, u, init_values.gamete_scheme, opti_G[t%L], s, dipl);

			old_distr = results.new_distr;

			// Test if the population fluctuation is stable from one cycle to the other
			if ((t%L == 0) & (L>0)){
	
				// Calculate the distribution distance between the two cycles;
				double distance = 0;
				for (int i = 0; i<ng; i++){
					distance += pow(old_distr[i] - data.all_distr[0][i], 2);
				}
				distance = sqrt(distance);

				// Check if the distance is smaller than 10^-12 in which case we stop the simulation
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

	// WITH STOCHASTICITY (meaning d>0)

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
			results = new_distributions(init_values.gamete_values, init_values.rec_table, old_distr, u, init_values.gamete_scheme, opti_G[t], s, dipl);
			data.all_distr[t%(L*10)] = results.new_distr;
			data.mean_fitness[t%(L*10)] = results.mean_fitness;
			old_distr = results.new_distr; 
		}
	}

	return data;
}

////////////////
// STATISTICS //
////////////////

/**
 * @brief Get the statistics of the article of Burger(2002) from the simulation.
 * 
 * @param gamete_distr 2D vecotr containing the distribution of each gamete at each timestep
 * @param gamete_values vector with the value of each gamete
 * @param loci_values value of each positive loci
 * @param mean_fitness vector of the mean fitness at each timestep
 * @param n the number of loci
 * @return a struct results containing:
 * 		- the averaged value of the mean genotype in the population
 * 		- the averaged genetic variance
 * 		- the ratio btw the variance and Vmax (the maximum possible variance)
 * 		- the geometric mean of the fitness 
 */
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

////////////////////
// INITIALISATION //
////////////////////

/**
 * @brief Generates the positive alleles values (of the diallelic loci)
 * 
 * @param dipl the population is diploid (= 1) or haploid (= 0)
 * @return init, is a struture that contains, the loci values, the value of each gamete, 
 * the structure of each gamete (i.e. bianry number with 0 = null_allele, 1 = positive_allele),
 * and the recombination table, a 3D vector where we can read 
 * the probability to have gamete i from gamete j and gamete k
 */ 
init initialisation(bool dipl)
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
	if (dipl == 1){
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

////////////////////////
// JK_I_RECOMBINATION //
////////////////////////

/**
 * @brief Gives the probability to obtain the gamete i through recombination of gamete j and k
 * 
 * @param i string corresponding to gamete i
 * @param j string corresponding to gamete j
 * @param k string corresponding to gamete k
 * @param r vector of the rate of recombination
 * @return the probability to get i from the j,k pair 
 */
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

///////////////////
// PERMUT_GAMETE //
///////////////////

/**
 * @brief Creates a new gamete from two gamete and the positions where a recombination occur
 * 
 * @param j the string corresponding to gamete j
 * @param k the string corresponding to gamete k
 * @param rec the vector with the positions where a recombination happens
 * @return the string of the new gamete after the recombination 
 */
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


///////////////////////
// NEW DISTRIBUTIONS //
///////////////////////

/**
 * @brief from the current distribution of gametes and the new optimum genotype,
 * calculates the new distribution of gametes in the population
 * 
 * @param gamete_values the vector containing the value of each gamete
 * @param dtf_rec the 3d vector with for each gametes j and k, the probability to have i through recombinatino
 * @param gamete_distr the current distribution of each gametes
 * @param mut_rate the mutation rate of the population
 * @param gam_bin the structure in binary of each gamete
 * @param opti_Gt the current optimum genotype for this environment
 * @param strength the strength of selection
 * 
 * @return a struct with the new distribution of gametes the mean fitness and the genetic variance
 */

one_step new_distributions(vector<double> gamete_values, vector<vector<vector<double>>> dtf_rec,
					   vector<double> gamete_distr, double mut_rate, vector<string> gam_bin, double opti_Gt, double strength, bool dipl)
{
	int ng = pow(2, n);
	vector<double> p_star(ng, 0);

	double sum = 0;
		for (int i = 0; i<ng; i++){
			sum += gamete_distr[i];
		}

	double mean_fitness = 0;
	
	
	if (dipl == 1){ // DIPLOID
		
		// Calculate the mean fitness
		for (int j = 0; j < ng; j++)
		{
			for (int k = j; k < ng; k++)
			{

				mean_fitness += exp(-strength * pow(gamete_values[j] + gamete_values[k] - opti_Gt, 2)) * gamete_distr[j] * gamete_distr[k];
			}
		}
		// Calculate the new distribution of gamete i (whitout mutation)
		for (int i = 0; i < ng; i++)
		{
			for (int j = 0; j < ng; j++)
			{
				for (int k = j; k < ng; k++) // need to consider every pair of parent gamete to add recombination
				{	
					p_star[i] += exp(-strength * pow(gamete_values[j] + gamete_values[k] - opti_Gt, 2)) * gamete_distr[j] * gamete_distr[k] * dtf_rec[i][j][k] / mean_fitness;
				}
			}
		}

	} else { // HAPLOID
		// Calculate the mean fitness
		for (int i = 0; i < ng; i++) 
			{
				mean_fitness += exp(-strength * pow(gamete_values[i]- opti_Gt, 2)) * gamete_distr[i];
			}
		// Calculate the new distribution of gamete i (whitout mutation)
		for (int i = 0; i < ng; i++){
			p_star[i] += exp(-strength * pow(gamete_values[i] - opti_Gt, 2)) * gamete_distr[i] / mean_fitness;
		}
	}
	
	

	// Adding mutation to the new distributions
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

