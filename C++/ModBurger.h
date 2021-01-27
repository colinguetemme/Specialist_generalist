#include<iostream> 
#include<chrono> 
#include<random> 
#include<bitset>
using namespace std; 
  
double *rnorm(int n, double mean, double sd) 
{ 
    /* Create random engine with the help of seed */
    unsigned seed = chrono::steady_clock::now().time_since_epoch().count(); 
    default_random_engine e (seed); 
  
    /* declaring normal distribution object 'distN' and initializing its mean and standard deviation fields. */
  /* Mean and standard deviation are distribution parameters of Normal distribution. Here, we have used mean=5, and standard deviation=2. You can take mean and standard deviation as per your choice */
    normal_distribution<double> distN(mean ,sd); 
    
    double list_rnorm[n];
    for (int i=1; i<=n; i++)
    {
      list_rnorm[i-1] = distN(e);
      cout<<i<<". "<<list_rnorm[i-1]<<"\n"; 
    }   
  double *p_list_rnorm;
  p_list_rnorm = &list_rnorm[0];
  std::cout << std::bitset<8>("11011011").to_ulong();
  return p_list_rnorm;  
}

