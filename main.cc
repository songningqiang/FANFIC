#include <iostream>

#include "prototype.h"

int main()
{

	//setup
	int Max_samp = 5e6; //number of sampling	
	std::string ordering = "NO"; //normal ordering or inverted ordering
	std::string t23oct = "upper"; //choose t23 octant from upper, lower and max, this only affects DUNE and HyperK	
	//list of experiments to use, if there are more than 1 exp with t23-dcp chi2 table, the last chi2 will be used
	//otherwise, the order of experiments doesn't matter
	//it is suggested to use one of DUNE, NUFIT and HYPERK with/without JUNO since they come with chi2 tables
	//if NUFIT is served, the NUFIT 5.0 t23-dcp chi2 table will be used
	//a complete list is {"JUNO", "NUFIT"/"HK"/"DUNE"}
	//if leave the list emtpy, the default NUFIT 5.0 table will be used, assuming each parameter is gaussian
	std::vector<std::string> experiments {"JUNO", "DUNE"};
	std::string foutput = "test.txt"; //name of output file to save flavor compositions and chi2
	bool randominitialflavor = false; //if true, sample initial flavor composition randomly from flat prior, otherwise, use the intial flavors below
	std::vector<double> initialflavor {1./3., 2./3., 0}; //initial flavor composition at the source, must sum up to 1	

	std::cout << Max_samp << " samples will be generated, with " << ordering << " ordering and " << t23oct << " octant." << std::endl;
	std::cout << "Output will be written in the file " << foutput << ", in the format of alpha_e, alpha_mu, alpha_tau, chi2." << std::endl; 
	oscillationparams osc(ordering); //this will initialize the oscillation parameters with NUFIT 5.0 table with superK by default

	//for JUNO
	JUNO JUNOEXP;

	//for DUNE
	std::string fname = "data/deltacp_theta23";
	if (t23oct == "upper") fname += "_";
	else if (t23oct == "lower") fname += "_woct_";
	else if (t23oct == "max") fname += "_max_";
	else {std::cout << "Wrong octant." << std::endl; exit(1);}		
	if (ordering == "NO") fname += "NH";
	else if (ordering == "IO") fname += "IH";
	else {std::cout << "Wrong mass ordering." << std::endl; exit(1);}
	fname += "_cursol.dat";
	DUNE DUNEEXP(fname);


	//for NUFIT 5.0 t23-dcp chi2 table
	if (ordering == "NO") fname = "data/nufit_5.0_skyes_no_s23sq_dcp.dat";
	else if (ordering == "IO") fname = "data/nufit_5.0_skyes_io_s23sq_dcp.dat";
	else {std::cout << "Wrong mass ordering." << std::endl; exit(1);}
	NUFIT NF(fname);

	//for HyperK t23-dcp chi2 table
	fname = "data/deltacp_theta23";
	if (t23oct == "upper") fname += "_";
	else if (t23oct == "lower") fname += "_woct_";
	else if (t23oct == "max") fname += "_max_";
	else {std::cout << "Wrong octant." << std::endl; exit(1);}		
	if (ordering == "NO") fname += "NH";
	else if (ordering == "IO") fname += "IH";
	else {std::cout << "Wrong mass ordering." << std::endl; exit(1);}
	fname += "_HK.dat";
	HYPERK HKEXP(fname);

	//list of experiments to use, if there are more than 1 exp with t23-dcp chi2 table, the last chi2 will be used
	//otherwise, the order of experiments doesn't matter
	//it is suggested to use one of DUNE, NUFIT and HYPERK with/without JUNO since they come with chi2 tables
	//if NUFIT is served, the NUFIT 5.0 t23-dcp chi2 table will be used
	//a complete list is {&JUNOEXP, &NF/&HKEXP/&DUNEEXP}
	//if leave the list emtpy, the default NUFIT 5.0 table will be used, assuming each parameter is gaussian
	//experimentlist explist = {&JUNOEXP, &NF};
	experimentlist explist;
	for (int i = 0; i < experiments.size(); ++i)
	{
		if (experiments[i] == "JUNO") {explist.push_back(&JUNOEXP); std::cout << "JUNO is included." << std::endl;}
		if (experiments[i] == "NUFIT") {explist.push_back(&NF); std::cout << "NUFIT 5.0 chi2 is included." << std::endl;}
		if (experiments[i] == "DUNE") {explist.push_back(&DUNEEXP); std::cout << "DUNE chi2 is included." << std::endl;}
		if (experiments[i] == "HK") {explist.push_back(&HKEXP); std::cout << "HYPERK chi2 is included." << std::endl;}
	}	
	likelihood oscchi2(explist, osc); //this will initialize the total chi2 module
	flavorregion flav; //this will initialize the flavor oscillation module	

	prior oscprior;
	FILE *fp;
 	fp = fopen(foutput.c_str(), "w");
 	fprintf(fp, "#nu_e nu_mu nu_tau chi2\n");
 	std:: cout << "Begin sampling..." << std::endl;
 	//start Monte Carlo
	for (int i = 0; i < Max_samp; ++i)
	{
		
		//draw random sin^2(theta_12) within 5sigma range
		double p_12 = oscprior.flatPrior(oscprior.rand01(),osc.get12best()-5.*osc.get12minus(), osc.get12best()+5.*osc.get12plus());
		double p_13 = oscprior.flatPrior(oscprior.rand01(),osc.get13best()-5.*osc.get13minus(), osc.get13best()+5.*osc.get13plus());
		//if t23-dcp chi2 data is available, use the data range there, dcp is in the unit of degrees
		double p_23, p_dcp;
		if (oscchi2.getchi2file() != "")
		{
			std::vector<double> v = oscchi2.gett23sqdata();
			p_23 = oscprior.flatPrior(oscprior.rand01(),*std::min_element(v.begin(), v.end()), 
				*std::max_element(v.begin(), v.end()));
			v = oscchi2.getdcpdata();
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
		if (randominitialflavor == true)
		{
			//draw random intial flavor composition
	      	comp_i[0] = oscprior.rand01();
	      	comp_i[1] = oscprior.flatPrior(oscprior.rand01(),0.,1.-comp_i[0]);
	      	comp_i[2] = 1.-comp_i[0]-comp_i[1];
	    }
	    else comp_i = initialflavor;
      	std::vector<double> comp_f = flav.evolveflavor(comp_i, oscp); //this will compute the flavor composition at the earth
      	//use t23-dcp chi2 to calculate the total chi2
      	double chi = 0.;
      	if (oscchi2.getchi2file() != "") chi = oscchi2.chisqfromdata(oscp);
      	else chi = oscchi2.chisq(oscp);
      	//save to file, keep only those within 5sigma for 2 dof
		if(chi < 30.) 
            fprintf(fp, "%1.4f %1.4f %1.4f %e\n", comp_f[0], comp_f[1], comp_f[2], chi);

	}
	fclose(fp);
	std:: cout << "Done sampling, exit." << std::endl;

}