#include "ModBurger.h"

int main(){

	string filename = "results/n";
	filename.append(to_string(n));
	filename.append("s");
	filename.append(to_string(int(s)));
	filename.append(".csv");
	cout << filename << endl;
	// + "n" + to_string(n) + "s" + to_string(s)
	//		+ "L" + to_string(L) + "d" + to_string(d) + "A" + to_string(A) + ".csv";

	std::ofstream outfile (filename);
	int tmax = 3000; //!// should be 300000 if d=0 or 50000 if d>0
	vector<double> opti_G = environment(tmax, L, d, A);
	
	int num_sim = 1;
	for (int e = 0; e<num_sim; e++){
		data = one_simul(0.5, 24, 0, 5, 0, 4, opti_G);
	}

	for (int i = 0; i<data.all_mean_fitness.size(); i++){
		cout << i << " : " << data.all_mean_fitness[i] << endl;
		cout << i << " : " << data.all_distr[i][1] << endl;
		cout << i << " : " << data.all_var_genetics[i] << endl;
	} 

	// creates the output files with the results

	outfile.close();
	return 0;
}

//------------------------

// function to do one simulation

many_steps one_simul(double A, int L, double d, double s, double u, int n, vector<double> opti_G){
	
	init init_values = initialisation();

	// Initialze the variable which will stock the distributions
	// of the gametes through the time
	one_step results;
	many_steps data;

	vector<vector<double>> all_distr(L, vector<double>(ng, 0));
	vector<double> mean_fit(L, 0);
	vector<double> var_gen(L, 0);

	// initialize the distribution
	data.all_distr = all_distr;
	data.all_mean_fitness = mean_fit;
	data.all_var_genetics = var_gen;
	vector<double> old_distr = init_values.gamete_distr;

	// If there is no stochasticity (meaning d==0)
	if (d==0){
		// Does the simulation
		for (int t = 0; t<opti_G.size(); t++){
			results = new_distributions(init_values.gamete_values, init_values.rec_table, old_distr, 0.00002, init_values.gamete_scheme, opti_G[t], s);

			old_distr = results.new_distr;

			// test whenever we need to check the stability
			if ((t%L == 0) & (L>0)){
	
				// Calculate the distance;
				double distance = 0;
				for (int i = 0; i<ng; i++){
					distance += pow(old_distr[i] - data.all_distr[0][i], 2);
				}
				distance = sqrt(distance);
				cout << t << " : " << distance << endl;

				// Check if the distance is higher than 10^-12

				if (distance < pow(10, -12)){
					return data;
				}
			}

			double mean_g = 0;
			double var_g = 0;
			for (int i = 0; i<ng; i++){
				mean_g += init_values.gamete_values[i]*results.new_distr[i];
			}
			for (int i = 0; i<ng; i++){
				var_g += pow((init_values.gamete_values[i]*results.new_distr[i])-mean_g, 2);
			}

			data.all_distr[t%L] = results.new_distr;
			data.all_mean_fitness[t%L] = results.mean_fitness;
			data.all_var_genetics[t%L] = var_g;

			double sum = 0;
			for (int i = 0; i<ng; i++){
				sum += data.all_distr[t%L][i];
			}
			
			// // outfile << t << "\t";
			// // outfile << opti_G[t] << "\t";
			// // for (int z = 0; z<ng; z++){
			// // 	outfile << all_distr[t+1][z] << "\t";
			// // }
			// // outfile << "\n";	
		}
	}

	// if there is stochasticity
	if (d>0){
		for (int t = 0; t<opti_G.size(); t++){
			results = new_distributions(init_values.gamete_values, init_values.rec_table, data.all_distr[t], 0.00002, init_values.gamete_scheme, opti_G[t], s);
			data.all_distr[t+1] = results.new_distr;
			data.all_mean_fitness[t] = results.mean_fitness;

			double mean_g = 0;
			double var_g = 0;
			for (int i = 0; i<ng; i++){
				mean_g += init_values.gamete_values[i]*results.new_distr[i];
			}
			for (int i = 0; i<ng; i++){
				var_g += pow((init_values.gamete_values[i]*results.new_distr[i])-mean_g, 2);
			}
			cout << var_g << endl;
			data.all_var_genetics[t] = var_g;
			// // outfile << t << "\t";
			// // outfile << opti_G[t] << "\t";
			// // for (int z = 0; z<ng; z++){
			// // 	outfile << all_distr[t+1][z] << "\t";
			// // }
			// // outfile << "\n";	
		}
	}

	return data;

}

// environment will define the optimum genotype value through time

vector<double> environment(int tmax, int period, double d, double L)
{
	opti_G.clear();
	double var = d * L;
	std::normal_distribution<> norm(0, var);
	for (int t = 0; t < tmax; t++)
	{
		opti_G.push_back(norm(rdgen) + 0.5 + L * sin(2 * M_PI * t / period));

	}
	return opti_G;
}

// initialisation will generate the phenotypic values of each loci,
// vector of recombination probabilities and the values of recombination
// probabilities from gamete j/k to gamete i.

init initialisation(void)
{	
	init init_values;
	int s, size, g;
	string bin;
	double sum = 0.0;
	//vector<double> alpha(n);
	alpha[0] = unif(rdgen);
	sum += alpha[0];

	for (int i = 0; i < n - 1; i++)
	{
		alpha[i + 1] = unif(rdgen); //initialise alleles
		sum += alpha[i + 1];
		r.push_back(0.0);
		r[i] = unif(rdgen) / 2.0; //initialise recombination rate
	}
	// cout << "Allelic values: " << endl;
	for (int i = 0; i < n; i++)
	{
		alpha[i] = 0.5 * alpha[i] / sum; //scale alleles

		// cout << alpha[i] << " ";
	}
	// cout << endl;

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

	init_values.gamete_values = genotypes;
	init_values.gamete_distr = distr;
	init_values.gamete_scheme = gametes;
	init_values.rec_table = rec_table;
	return init_values;
}

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
	//cout << R << endl;
	return R;
}

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
	//cout << opti_Gt << endl
	int ng = pow(2, n);
	vector<double> p_star(ng, 0);
	double mean_fitness = 0;
	double sum = 0;
		for (int i = 0; i<ng; i++){
			sum += gamete_distr[i];
		}

		
	// cout << "s = " << strength << endl;
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

