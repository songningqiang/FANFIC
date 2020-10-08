#include <iostream>

#include "prototype.h"

int main()
{

	//setup
	int Max_samp = 1e7; //number of sampling
	int year = 2020; //can only choose 2020 (current) or 2040 (future)

	bool normalized = true; //if true, will normalize the sum of final flavor compositions to be 1
	bool randominitialflavor = true; //if true, sample initial flavor/mass composition randomly from flat prior, otherwise, use the intial flavors/masses below
	std::vector<double> initialflavor {1./3., 2./3., 0.}; //initial flavor composition at the source, must sum up to 1
	std::string foutput = "test.txt"; //name of output file to save flavor compositions and chi2

	std::cout << "Output will be written in the file " << foutput << ", in the format of alpha_e, alpha_mu, alpha_tau, chi2." << std::endl; 
	

	nonunitflavorregion flav(year); //this will initialize the flavor oscillation module
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
		if (randominitialflavor == true)
		{
			//draw random intial flavor composition
	      	comp_i[0] = oscprior.rand01();
	      	comp_i[1] = oscprior.flatPrior(oscprior.rand01(),0.,1.-comp_i[0]);
	      	comp_i[2] = 1.-comp_i[0]-comp_i[1];
	    }
	    else comp_i = initialflavor;
      	std::vector<double> comp_f(3, 0.);
      	comp_f= flav.evolvefromflavor(comp_i, Usqinput, normalized); //this will compute the flavor composition at the earth from a combination of mass states
      	//use t23-dcp chi2 to calculate the total chi2
      	double chi = flav.chisq(Usqinput);
      	//save to file, keep only those within 5sigma for 2 dof
		if(chi < 30.) 
            fprintf(fp, "%1.4f %1.4f %1.4f %e\n", comp_f[0], comp_f[1], comp_f[2], chi);

	}
	fclose(fp);
	std:: cout << "Done sampling, exit." << std::endl;	


}