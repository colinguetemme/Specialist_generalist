#include "ModBurger.h"

int main() {



	int ngg = pow(2.0, nn);

	int** t;

	t = new int* [ngg];
	for (int z = 0; z < ngg; z++) {
		t[z] = new int[nn];
		bitset<nn> allele(z);
		for (int zz = 0; zz < nn; zz++) {
			t[z][zz] = allele[zz] * (zz+1);
			cout << t[z][zz] << " ";
		}
		cout << endl;
	}

	string p = "3";
	string g = "6";
	string gg = p + g;
	cout << gg;


	//initialisation();
	//rec_test();
	int x;
	cin >> x;

	return 0;
}

//------------------------
void initialisation(void) {
	int s, size, g;
	double sum = 0.0;
	string bin;

	alpha[0] = unif(rdgen);
	sum += alpha[0];

	for (int i = 0; i < n-1; i++) {
		alpha[i+1] = unif(rdgen); //initialise alleles
		sum += alpha[i+1];

		r[i] = unif(rdgen) / 2.0; //initialise recombination rate
	}
	cout << "Allelic values: " << endl;
	for (int i = 0; i < n; i++) {
		alpha[i] = 0.5 * alpha[i] / sum; //scale alleles

		cout << alpha[i] << " ";
	}
	cout << endl;

	cout << "Gametes: " << endl;

	//initialise gametes
	for (int i = 0; i < ng; i++) {
		
        bitset<n> alleles(i);
		gametes.push_back(alleles.to_string<char,std::string::traits_type,std::string::allocator_type>());
		sum = 0.0;

		for (int j = 0; j < n; j++) {			
			if (((int)gametes[i][j])-48 == 1) sum += alpha[j];
		}
		genotypes.push_back(sum);
		cout << sum << "\t";
		
	}

	//initiliase recombination table
	rec_table = new double **[ng];
	for (int j = 0; j < ng; j++) {
		rec_table[j] = new double *[ng];
		for (int k = 0; k < ng; k++) {
			rec_table[j][k] = new double[ng];
			for (int i = 0; i < ng; i++) {
				rec_table[j][k][i] = jk_i_recombination(gametes[i], gametes[j], gametes[k], r);
			}
		}

	}

}



double jk_i_recombination(string i, string j, string k, vector<double> rec) {

	double R = 0.0;
	double prod = 1.0;
	int ngg = pow(2.0, nn);
	int sum;
	int** t;

	vector<int> rec_vec;
	vector<string> jk;

	t = new int *[ngg];
	for (int z = 0; z < ngg; z++) {
		t[z] = new int[nn];
		bitset<nn> allele(z);
		sum = 0;
		for (int zz = 0; zz < nn; zz++) {
			t[z][zz] = allele[zz] * (zz + 1);
			if (t[z][zz] > 0) rec_vec.push_back(t[z][zz]);
		}

		if (z == 0) jk = permut_gamete(j, k, rec_vec);
		else {
			jk.push_back(j);
			jk.push_back(k);
			for (int zz = 1; zz < nn; zz++) rec_vec.push_back(-zz);
		}

		for (int a = 0; a < 2; a++) {
			if (jk[a] == i) {
				int rv = 0;
				prod = 1.0;
				for (int zz = 0; zz < (int)r.size(); zz++) {
					if (zz == rec_vec[rv] - 1) {
						prod *= r[zz];
						rv++;
					}
					else prod *= 1.0 - r[zz];
				}
				R += 0.5 * prod;
			}
		}


	}

	return R;
}

vector<string> permut_gamete(string j, string k, vector<int> rec) {
	int rec_spot;

	string j_bin_rec, k_bin_rec;
	vector<string> jk;

	rec_spot = rec[0];

	for (int z = 0; z < rec_spot; z++) j_bin_rec += j[z];
	for (int z = rec_spot; z < (int)j.size(); z++) j_bin_rec += k[z];

	for (int z = 0; z < rec_spot; z++) k_bin_rec += k[z];
	for (int z = rec_spot; z < (int)k.size(); z++) k_bin_rec +=j[z];

	rec.erase(rec.begin());

	if ((int)rec.size() > 1) permut_gamete(j_bin_rec, k_bin_rec, rec);

	jk.push_back(j_bin_rec);
	jk.push_back(k_bin_rec);

	return jk;
}













double recombination(string j, string k, string i) {
	double R = 1.0;
	int gam; //j = 0, k = 1

	if (i[0] == j[0] || i[0] == k[0]) {

		if ((i[0] == j[0] && i[0] == k[0]) || (i[0] == j[0] && i[0] != k[0])) gam = 0;
		else gam = 1;

		for (int z = 1; z < n; z++) {
			if (i[z] != j[z] && i[z] != k[z]) {
				R = 0.0;
				break;
			}
			if (i[z] == j[z] && i[z] == k[z]) {
				double rdn = unif(rdgen);
				if (rdn < r[z]) {
					if (gam == 0) gam = 1;
					else gam = 0;
				}
			}
			else {
				if (gam == 0) { //j is the starting gamete				
					if (i[z] != j[z] && i[z] == k[z]) {
						R *= r[z];
						gam = 1;
					}
					if (i[z] == j[z] && i[z] != k[z]) R *= (1.0 - r[z]);
				}
				else { //k is the starting gamete
					if (i[z] == j[z] && i[z] != k[z]) {
						R *= r[z];
						gam = 0;
					}
					if (i[z] == j[z] && i[z] != k[z]) R *= (1.0 - r[z]);
				}
			}
		}

	}
	else R = 0.0;
	return R;
}




	

/*
rec_table recombination_table(int i){
	for (int j = 0; j < n; j++){
		for (int k = j; k < n; k++){
			for (int l = 0; l < (n-1); l++){
				
			}
		}
	}
}
*/
