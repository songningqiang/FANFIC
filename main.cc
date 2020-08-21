#include <iostream>

#include "prototype.h"

int main()
{

	//setup
	std::string ordering = "NO";
	int Max_samp = 5e6; //number of sampling
	std::string foutput = "test.txt"; //name of output file to save flavor compositions and chi2
	bool useJUNO = true;
	bool useDUNE = true;
	bool useNUFIT = false; //use SM t23-dcp chi2 or not
	bool useHK = false; //use hyperK t23-dcp chi2 or not
	std::string t23oct = "upper"; //choose t23 octant from upper, lower and max
	oscillationparams osc(ordering); //this will initialize the oscillation parameters with NUFIT 5.0 table with superK by default
	JUNO JUNOEXP(osc);
	DUNE *DUNEp;
	HYPERK *HKp;
	NUFIT *NFp;
	if (useDUNE)
	{
		std::string fname = "data/deltacp_theta23";
		if (t23oct == "upper") fname += "_";
		else if (t23oct == "lower") fname += "_woct_";
		else if (t23oct == "max") fname += "_max_";
		else {std::cout << "Wrong octant." << std::endl; exit(1);}		
		if (ordering == "NO") fname += "NH";
		else if (ordering == "IO") fname += "IH";
		else {std::cout << "Wrong mass ordering." << std::endl; exit(1);}
		fname += "_cursol.dat";
		DUNE *DUNEEXP = new DUNE(osc, fname);
		DUNEp = DUNEEXP;	
	}

	if (useNUFIT)
	{
	
		if (ordering == "NO") {NUFIT *NF = new NUFIT(osc, "nufit_5.0_skyes_no_s23sq_dcp.dat"); NFp = NF;}
		else if (ordering == "IO") {NUFIT *NF = new NUFIT(osc, "nufit_5.0_skyes_io_s23sq_dcp.dat"); NFp = NF;}
		else {std::cout << "Wrong mass ordering." << std::endl; exit(1);}
	}

	if (useHK)
	{
		std::string fname = "data/deltacp_theta23";
		if (t23oct == "upper") fname += "_";
		else if (t23oct == "lower") fname += "_woct_";
		else if (t23oct == "max") fname += "_max_";
		else {std::cout << "Wrong octant." << std::endl; exit(1);}		
		if (ordering == "NO") fname += "NH";
		else if (ordering == "IO") fname += "IH";
		else {std::cout << "Wrong mass ordering." << std::endl; exit(1);}
		fname += "_HK.dat";
		HYPERK *HKEXP = new HYPERK(osc, fname);
		HKp = HKEXP;	
	}

	//list of experiments to use, if there are more than 1 exp with t23-dcp chi2, the last chi2 will be used
	experimentlist explist = {&JUNOEXP, DUNEp};
	//experimentlist explist = {};
	oscillationprob oscprob(explist, osc); //this will initialize the total chi2 module
	flavorregion flav; //this will initialize the flavor oscillation module	

	prior oscprior;
	FILE *fp;
 	fp = fopen(foutput.c_str(), "w");
	for (int i = 0; i < Max_samp; ++i)
	{
		
		double p_12 = oscprior.flatPrior(oscprior.rand01(),osc.get12best()-5.*osc.get12minus(), osc.get12best()+5.*osc.get12plus());
		double p_13 = oscprior.flatPrior(oscprior.rand01(),osc.get13best()-5.*osc.get13minus(), osc.get13best()+5.*osc.get13plus());
		//if t23-dcp chi2 data is available, use the data range there, dcp is in the unit of degrees
		double p_23, p_dcp;
		if (oscprob.getchi2file() != "")
		{
			std::vector<double> v = oscprob.gett23sqdata();
			p_23 = oscprior.flatPrior(oscprior.rand01(),*std::min_element(v.begin(), v.end()), 
				*std::max_element(v.begin(), v.end()));
			v = oscprob.getdcpdata();
			p_dcp = oscprior.flatPrior(oscprior.rand01(),*std::min_element(v.begin(), v.end()), 
				*std::max_element(v.begin(), v.end()));			
		}
		else
		{
			p_23 = oscprior.flatPrior(oscprior.rand01(),osc.get23best()-5.*osc.get23minus(), osc.get23best()+5.*osc.get23plus());
			p_dcp = oscprior.flatPrior(oscprior.rand01(), 0., 540.);
		}
		//oscillation parameters in the order of \sin^2\theta_{12}, \sin^2\theta_{13}, \sin^2\theta_{23}, dcp
		//dcp is in the unit of degrees						
		std::vector<double> oscp{p_12, p_13, p_23, p_dcp};
		//initial flavor composition in the order of \nu_\mu, \nu_e, \nu_\tau
		std::vector<double> comp_i(3,0.);
      	comp_i[0] = oscprior.rand01();
      	comp_i[1] = oscprior.flatPrior(oscprior.rand01(),0.,1.-comp_i[0]);
      	comp_i[2] = 1.-comp_i[0]-comp_i[1];
      	std::vector<double> comp_f = flav.evolveflavor(comp_i, oscp); //this will compute the flavor composition at the earth
      	//use t23-dcp chi2 to calculate the total chi2
      	double chi = 0.;
      	if (oscprob.getchi2file() != "") chi = oscprob.chisqfromdata(oscp);
      	else chi = oscprob.chisq(oscp);
      	//save to file
		if(chi < 30.) 
            fprintf(fp, "%1.4f %1.4f %1.4f %e\n", comp_f[0], comp_f[1], comp_f[2], chi);

	}
	fclose(fp);

}