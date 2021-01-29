#include "ModBurger.h"

int main() {

	initialisation();
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
		
        std::bitset<n> alleles(i);
		gametes.push_back(alleles.to_string<char,std::string::traits_type,std::string::allocator_type>());
		sum = 0.0;

		for (int j = 0; j < n; j++) {			
			if (((int)gametes[i][j])-48 == 1) sum += alpha[j];
		}
		genotypes.push_back(sum);
		cout << sum << "\t";
		
	}

}
