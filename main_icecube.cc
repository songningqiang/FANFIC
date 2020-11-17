#include <iostream>

#include "prototype.h"

/************************************************************************************************************/
/*get the posterior of x in the source flavor composition (x:1-x:0) with both oscillation and Icecube data*/

int main()
{

	//setup
	int Max_samp = 5e7; //number of sampling
	double nufitversion = 5.0; //choose nufit version from {1.0,1.1,1.2,1.3,2.0,2.1,2.2,3.0,3.1,3.2,4.0,4.1,5.0}
	std::string	year = "2015"; //IceCube year, choose 2015 (NUFIT), 2020 (NUFIT), 2028 (NUFIT+JUNO) or 2040 (JUNO+DUNE+HK)
	std::string ordering = "IO"; //normal ordering or inverted ordering
	std::string t23oct = "upper"; //choose t23 octant from upper, lower and max, this only affects DUNE and HyperK	
	//list of experiments to use, if there are more than 1 exp with t23-dcp chi2 table, the sum of chi2 will be used
	//in that case, make sure the chi2 tables match each other, otherwise an error will be raised
	//the order of experiments in the list doesn't matter
	//it is suggested to use one of DUNE, NUFIT and HYPERK with/without JUNO since they come with chi2 tables
	//if NUFIT is served, the NUFIT vx.x t23-dcp chi2 table will be used if version>=2.0
	//a complete list is {"JUNO", "NUFIT", "HK", "DUNE"}
	//if leave the list emtpy, no chi2 table will be used, each parameter is assumed to be gaussian with centers and errors specified by the Nufit version
	//std::vector<std::string> experiments {"JUNO", "DUNE", "HK"};
	std::vector<std::string> experiments {"NUFIT"};
	//Use user-specified dcp value or not, can be true only if HK or DUNE is used and upper octant is specified
	//if false, the default Nufit5.0 bestfit will be implemented
	//if false, 1.1pi is used for NO and 1.5pi is used for IO
	bool customizedcp = false;
	//True value of dcp in units of pi, this option is only available for HK and DUNE and upper octant
	//For NO, choose between 0, 0.5 and 1, 1 corresponds to the Nufit5.0 bestfit of 1.1pi
	//For IO, choose between 0, 1 and 1.5, 1.5 corresponds to the Nufit5.0 bestfit of 1.5pi
	double dcpbest = 1.;
	std::string foutput = "test.txt"; //name of output file to save flavor compositions and chi2	
	//std::string foutput = "output_testsource/output_NUFIT5.0_IO_IceCube2015.txt";



	oscillationparams osc(ordering, nufitversion); //this will initialize the oscillation parameters with NUFIT 5.0 table with superK by default
	std::cout << "Using Nufit version " << osc.getversion() << std::endl;
	std::cout << "Icecube year: " << year << std::endl;
	std::cout << Max_samp << " samples will be generated, with " << ordering << " ordering and " << t23oct << " octant." << std::endl;
	if (customizedcp == true) std::cout << "Use customized delta_CP = " << dcpbest << "*pi." << std::endl;
	else std::cout << "Use delta_CP from the bestfit of Nufit" << osc.getversion() << "." << std::endl;
	std::cout << "Output will be written in the file " << foutput << ", in the format of alphe_e (source), alpha_e, alpha_mu, alpha_tau, chi2." << std::endl; 


	//for JUNO
	JUNO JUNOEXP;

	//for DUNE
	std::string fname = "data/deltacp_theta23sq_Chi2/deltacp_theta23";
	if (t23oct == "upper") fname += "_";
	else if (t23oct == "lower") fname += "_woct_";
	else if (t23oct == "max") fname += "_max_";
	else {std::cout << "Wrong octant." << std::endl; exit(1);}		
	if (ordering == "NO") fname += "NH";
	else if (ordering == "IO") fname += "IH";
	else {std::cout << "Wrong mass ordering." << std::endl; exit(1);}	
	fname += "_cursol";
	if (customizedcp == true)
	{
		if (ordering == "NO")
		{
			if (dcpbest == 0.) fname += "_delta0";
			else if (dcpbest == 0.5) fname += "_delta0.5";
			else if (dcpbest == 1.) fname += "";
			else {std::cout << "Wrong dcp specified, must choose in between 0, 0.5 and 1 for NO." << std::endl; exit(1);}
		}
		if (ordering == "IO")
		{
			if (dcpbest == 0.) fname += "_delta0";
			else if (dcpbest == 1.) fname += "_delta1.0";
			else if (dcpbest == 1.5) fname += "";
			else {std::cout << "Wrong dcp specified, must choose in between 0, 1 and 1.5 for IO." << std::endl; exit(1);}
		}
	}			
	fname += ".dat";
	DUNE DUNEEXP(fname);


	//for NUFIT vx.x t23-dcp chi2 table
	fname = "";
	if (nufitversion >= 2.0)
	{
		fname += "data/deltacp_theta23sq_Chi2/v"+std::to_string(int(nufitversion*10))+".release-data-"+ordering+"_s23sq_dcp.dat";
		if (std::find(experiments.begin(), experiments.end(), "NUFIT") != experiments.end())
			std::cout << "Nufit chi2 table " << fname << " is being used." << std::endl;
	}
	NUFIT NF(fname);

	//for HyperK t23-dcp chi2 table
	fname = "data/deltacp_theta23sq_Chi2/deltacp_theta23";
	if (t23oct == "upper") fname += "_";
	else if (t23oct == "lower") fname += "_woct_";
	else if (t23oct == "max") fname += "_max_";
	else {std::cout << "Wrong octant." << std::endl; exit(1);}		
	if (ordering == "NO") fname += "NH";
	else if (ordering == "IO") fname += "IH";
	else {std::cout << "Wrong mass ordering." << std::endl; exit(1);}
	fname += "_HK";
	if (customizedcp == true)
	{
		if (ordering == "NO")
		{
			if (dcpbest == 0.) fname += "_delta0";
			else if (dcpbest == 0.5) fname += "_delta0.5";
			else if (dcpbest == 1.) fname += "";
			else {std::cout << "Wrong dcp specified, must choose in between 0, 0.5 and 1 for NO." << std::endl; exit(1);}
		}
		if (ordering == "IO")
		{
			if (dcpbest == 0.) fname += "_delta0";
			else if (dcpbest == 1.) fname += "_delta1.0";
			else if (dcpbest == 1.5) fname += "";
			else {std::cout << "Wrong dcp specified, must choose in between 0, 1 and 1.5 for IO." << std::endl; exit(1);}
		}
	}		
	fname += ".dat";
	HYPERK HKEXP(fname);


	//list of experiments to use, if there are more than 1 exp with t23-dcp chi2 table, the sum of chi2 will be used
	//in that case, make sure the chi2 tables match each other, otherwise an error will be raised
	//the order of experiments in the list doesn't matter
	//it is suggested to use one of DUNE, NUFIT and HYPERK with/without JUNO since they come with chi2 tables
	//if NUFIT is served, the NUFIT vx.x t23-dcp chi2 table will be used if version>=2.0
	//a complete list is {&JUNOEXP, &NF, &HKEXP, &DUNEEXP}
	//if leave the list emtpy, no chi2 table will be used, each parameter is assumed to be gaussian with centers and errors specified by the Nufit version
	//experimentlist explist = {&JUNOEXP, &NF};
	experimentlist explist;
	for (int i = 0; i < experiments.size(); ++i)
	{
		if (experiments[i] == "JUNO") {explist.push_back(&JUNOEXP); std::cout << "JUNO is included." << std::endl;}
		if (experiments[i] == "NUFIT") {explist.push_back(&NF); std::cout << "NUFIT" << osc.getversion() << " chi2 is included." << std::endl;}
		if (experiments[i] == "DUNE") {explist.push_back(&DUNEEXP); std::cout << "DUNE chi2 is included." << std::endl;}
		if (experiments[i] == "HK") {explist.push_back(&HKEXP); std::cout << "HYPERK chi2 is included." << std::endl;}
	}	
	likelihood oscchi2(explist, osc); //this will initialize the total chi2 module
	flavorregion flav; //this will initialize the flavor oscillation module

	likelihood_ice icechi2(year); ////this will initialize the IceCube chi2 module

	prior oscprior;
	FILE *fp;
 	fp = fopen(foutput.c_str(), "w");
 	fprintf(fp, "#x=f_e@source nu_e nu_mu nu_tau chi2\n");
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

			p_23 = oscprior.flatPrior(oscprior.rand01(), oscchi2.gett23sqdatamin(), oscchi2.gett23sqdatamax());
			p_dcp = oscprior.flatPrior(oscprior.rand01(), oscchi2.getdcpdatamin(), oscchi2.getdcpdatamax());			
		}
		else
		{
			p_23 = oscprior.flatPrior(oscprior.rand01(),osc.get23best()-5.*osc.get23minus(), osc.get23best()+5.*osc.get23plus());
			p_dcp = oscprior.flatPrior(oscprior.rand01(), 0., 540.);
		}
		

		//p_12 = 0.304; p_13 = 0.02219; p_23 = 0.573; p_dcp = 197.;
		//p_12 = 0.307; p_13 = 0.02206; p_23 = 0.538; p_dcp = 234.;
		//oscillation parameters in the order of \sin^2\theta_{12}, \sin^2\theta_{13}, \sin^2\theta_{23}, dcp
		//dcp is in the unit of degrees						
		std::vector<double> oscp{p_12, p_13, p_23, p_dcp};
		//initial flavor composition in the order of \nu_\mu, \nu_e, \nu_\tau
		std::vector<double> comp_i(3,0.);
      	comp_i[0] = oscprior.rand01();
      	comp_i[1] = 1.-comp_i[0];
      	comp_i[2] = 0.;

      	std::vector<double> comp_f(3, 0.);
      	comp_f= flav.evolvefromflavor(comp_i, oscp); //this will compute the flavor composition at the earth
      	
      	//use t23-dcp chi2 to calculate the total chi2
      	double chi = 0.;
      	if (oscchi2.getchi2file() != "") chi = oscchi2.chisqfromdata(oscp);
      	else chi = oscchi2.chisq(oscp);
      	chi += icechi2.chisqice(comp_f);
      	//save to file, keep only those within 5sigma for 2 dof
		if(chi < 30.) 
            fprintf(fp, "%1.6f %1.4f %1.4f %1.4f %e\n", comp_i[0], comp_f[0], comp_f[1], comp_f[2], chi);

	}
	fclose(fp);
	std:: cout << "Done sampling, exit." << std::endl;

}