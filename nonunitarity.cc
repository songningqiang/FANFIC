/***********************************************************************************************************/
//This is a Monte Carlo code to find the flavor compositions at Earth from certain source flavor compositions
//assuming nonunitarity  or future oscillation experiments
//f_{\beta,\oplus} = \sum_{i,\alpha} 1/(N_{\alpha}*N_{\beta})|U_{\alpha i}|^2*|U_{\beta i}|^2*f_{\alpha,S}
//Mar 19, 2021, Ningqiang Song
/***********************************************************************************************************/

#include <iostream>
//for parallelization
#include <omp.h>
#include "prototype.h"

int main()
{

	//setup
	int Max_samp = 5e7; //number of sampling
	int year = 2040; //can only choose 2020 (current) or 2040 (future)
	std::string oscoption = "submatrix"; //non-unitarity scheme, choose submatrix or agnostic, see arXiv:2008.01088 for details

	bool normalized = true; //if true, will normalize the sum of final flavor compositions to be 1
	bool randominitialflavor = true; //if true, sample initial flavor composition randomly from flat prior, otherwise, use the intial flavor below
	std::vector<double> initialflavor {0., 1., 0.}; //initial flavor composition at the source, must sum up to 1
	std::string foutput = "test.txt"; //name of output file to save flavor compositions and chi2

	std::cout << "Output will be written in the file " << foutput << ", in the format of alpha_e, alpha_mu, alpha_tau, chi2." << std::endl; 
	

	nonunitflavorregion flav(year, oscoption); //this will initialize the flavor oscillation module

	//print the best-fit of Usq matrix
	//printmatrix(flav.getUsqbest(), 3, 3);

	std::vector<double> Usqmin, Usqmax;	
	for (auto x : flav.getUsqdata())
	{
		Usqmin.push_back(*std::min_element(x.Usqdata.begin(), x.Usqdata.end()));
		Usqmax.push_back(*std::max_element(x.Usqdata.begin(), x.Usqdata.end()));
	}

	prior oscprior;
	FILE *fp;
 	fp = fopen(foutput.c_str(), "w");
 	fprintf(fp, "#nu_e nu_mu nu_tau chi2\n");
 	std:: cout << "Begin sampling..." << std::endl;
 	std::vector<double> Usqinput (3*3, 0.);
 	//start Monte Carlo
 	//uncomment this line if parallelization with openmp is necessary
 	//#pragma omp parallel for  	
	for (int i = 0; i < Max_samp; ++i)
	{
				
		//draw random |U_{ij}|^2 within parameter range
		std::vector<matrixdata> v = flav.getUsqdata();
		for (int k = 0; k < v.size(); ++k)
		{
			Usqinput[k] = oscprior.flatPrior(oscprior.rand01(), Usqmin[k], Usqmax[k]);
		}
		//initial flavor composition in the order of \nu_\mu, \nu_e, \nu_\tau
		std::vector<double> comp_i(3,0.);
		if (randominitialflavor == true) comp_i = oscprior.randomInitialFlavor();
	    else comp_i = initialflavor;
      	std::vector<double> comp_f(3, 0.);
      	comp_f= flav.evolvefromflavor(comp_i, Usqinput, normalized); //this will compute the flavor composition at the earth from a combination of mass states
      	//use t23-dcp chi2 to calculate the total chi2
      	double chi = flav.chisq(Usqinput);
      	//save to file, keep only those within 5sigma for 2 dof
		if(chi < 20.) 
            fprintf(fp, "%1.4f %1.4f %1.4f %e\n", comp_f[0], comp_f[1], comp_f[2], chi);

	}
	fclose(fp);
	std:: cout << "Done sampling, exit." << std::endl;	


}