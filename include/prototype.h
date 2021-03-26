#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <functional>
#include <chrono>
#include <ctime>
#include <random>

//special functions, needed for neutrino decay
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/random.hpp>


//for linear interpolation
//credit: Ronaldo Carpio
#include "linterp.h"
//for Dirichlet random numbers
//credit: George Cantwell https://github.com/gcant/dirichlet-cpp
#include "dirichlet.h"


//remove duplicated elements
template <typename ForwardIterator>
ForwardIterator remove_duplicates( ForwardIterator first, ForwardIterator last )
{
    auto new_last = first;

    for ( auto current = first; current != last; ++current )
    {
        if ( std::find( first, new_last, *current ) == new_last )
        {
            if ( new_last != current ) *new_last = *current;
            ++new_last;
        }
    }

    return new_last;
}

//print an m * n matrix from a 1d vector
void printmatrix(std::vector<double> v, int m, int n);

class oscillationparams
{
	public:
		oscillationparams(const std::string ordering = "NO", const double = 5.0);
		~oscillationparams(){};

		void set12(double t12sqb, double t12sqsigmap = NAN, double t12sqsigmam = NAN);
		void set13(double t13sqb, double t13sqsigmap = NAN, double t13sqsigmam = NAN);
		void set23(double t23sqb, double t23sqsigmap = NAN, double t23sqsigmam = NAN);
		void setdm21(double dm21sqb, double dm21sqsigmap = NAN, double dm21sqsigmam = NAN);
		void setdm31(double dm31sqb, double dm31sqsigmap = NAN, double dm31sqsigmam = NAN);
		void setdcp(double dcpb, double dcpsigmap = NAN, double dcpsigmam = NAN);
		double get12best() {return t12sqbest;}
		double get12plus() {return t12sqsigmaplus;}
		double get12minus() {return t12sqsigmaminus;}
		double get13best() {return t13sqbest;}
		double get13plus() {return t13sqsigmaplus;}
		double get13minus() {return t13sqsigmaminus;}
		double get23best() {return t23sqbest;}
		double get23plus() {return t23sqsigmaplus;}
		double get23minus() {return t23sqsigmaminus;}
		double getdm21best() {return dm21sqbest;}
		double getdm21plus() {return dm21sqsigmaplus;}
		double getdm21minus() {return dm21sqsigmaminus;}
		double getdm31best() {return dm31sqbest;}
		double getdm31plus() {return dm31sqsigmaplus;}
		double getdm31minus() {return dm31sqsigmaminus;}
		double getdcpbest() {return dcpbest;}
		double getdcpplus() {return dcpsigmaplus;}
		double getdcpminus() {return dcpsigmaminus;}
		std::string getordering() {return ordering;}
		std::string getversion() {return std::to_string(version);}

	private:
		double t12sqbest, t13sqbest, t23sqbest, dm21sqbest, dm31sqbest, dcpbest;
		double t12sqsigmaplus, t13sqsigmaplus, t23sqsigmaplus, dm21sqsigmaplus, dm31sqsigmaplus, dcpsigmaplus;
		double t12sqsigmaminus, t13sqsigmaminus, t23sqsigmaminus, dm21sqsigmaminus, dm31sqsigmaminus, dcpsigmaminus;
		std::string ordering;
		double version;
};


class matrixdata
{
	public:
		std::vector<double> Usqdata, Usqchi2;	
};

//in case of non-unitarity, parameterize the mixing matrix directly
class nonunitflavorregion
{
	public:
		nonunitflavorregion(const int year = 2020, const std::string oscoption = "submatrix");
		~nonunitflavorregion(){};
		void readchi2(std::string chi2file);
		void fillchi2(const int year, const std::string oscoption);
		double chisq(std::vector<double> Usqinput);
		std::vector<double> evolvefromflavor(std::vector<double> comp_i, std::vector<double> Usqinput, bool normalized = false);
		std::vector<matrixdata> getUsqdata() {return Usq;}
		//2d array of pointers pointing to the interpolation function of each matrix element and its chi2
		std::vector<double> getUsqbest() {return Usqbest;}
		InterpMultilinear<1, double> *interp1d[3*3];		
	private:
		std::vector<matrixdata> Usq;
		std::vector<double> Usqbest;
		std::vector<double> comp_f {std::vector<double>(3,0.)};	
	
};


class flavorregion
{
	public:
		flavorregion(){};
		~flavorregion(){};
		double sq(double x);
		double V11s(double q12, double q13);
		double V12s(double q12, double q13);
		double V13s(double q13);
		double V21s(double q12, double q13, double q23, double dcp);
		double V22s(double q12, double q13, double q23, double dcp);
		double V23s(double q13, double q23);
		double V31s(double q12, double q13, double q23, double dcp);
		double V32s(double q12, double q13, double q23, double dcp);
		double V33s(double q13, double q23);
		std::vector<double> evolvefromflavor(std::vector<double> comp_i, std::vector<double> oscinput);
		std::vector<double> evolvebackfromflavor(std::vector<double> comp_f, std::vector<double> oscinput);
		std::vector<double> evolvefrommass(std::vector<double> comp_i, std::vector<double> oscinput);
		double *cholesky(double *A, int n);
		double *lowermatrixinversion(double *A, int n);
		double *symmetricmatrixinversion(double *A, int n);
		void show_matrix(double *A, int n);
	private:
		//std::vector<double> comp_f {std::vector<double>(3,0.)};	
};


class oscillationexperiment
{
	public:
		oscillationexperiment(const std::string chi2file = "");
		~oscillationexperiment(){};
		virtual void setosc(oscillationparams &osc){}; //set oscillation parameters in an experiment
		void readchi2(std::string chi2file);
		std::vector<double> gett23sqdata() {return t23sqdata;}
		std::vector<double> getdcpdata() {return dcpdata;}
		std::vector<double> getchi2data() {return chi2data;}
		std::string getchi2file() {return schi2file;}
	protected:
		std::vector<double> t23sqdata, dcpdata, chi2data;
		std::string schi2file = "";
};

typedef std::vector<oscillationexperiment*> experimentlist;

class likelihood
{
	public:
		likelihood(experimentlist explist, oscillationparams &oscparams);
		~likelihood(){};
		void oscexperiments(experimentlist explist, oscillationparams &osc);
		double sq(double x);	
		double chisq(std::vector<double> oscinput);
		double chisq1213(std::vector<double> oscinput);
		double chisq23dcp(std::vector<double> oscinput);
		double chisqfromdata(std::vector<double> oscinput);
		std::vector<double> gett23sqdata() {return t23sqdata;}
		std::vector<double> getdcpdata() {return dcpdata;}
		std::vector<double> getchi2data() {return chi2data;}
		double gett23sqdatamin() {return t23sqdatamin;}
		double gett23sqdatamax() {return t23sqdatamax;}
		double getdcpdatamin() {return dcpdatamin;}
		double getdcpdatamax() {return dcpdatamax;}
		std::string getchi2file() {return schi2file;}	
	private:
		std::vector<double> t23sqdata, dcpdata, chi2data;
		double t23sqdatamin = 0., t23sqdatamax = 0., dcpdatamin = 0., dcpdatamax = 0.;
		std::string schi2file = "";
		oscillationparams osc;
		InterpMultilinear<2, double> *interp2d;		
};


class likelihood_ice
{
	public:
		likelihood_ice(std::string year);
		~likelihood_ice(){};

		void readchi2(std::string chi2file);
		double chisqice(std::vector<double> flavinput);
		double chisqice_2015(std::vector<double> flavinput);
		std::vector<double> getfedata() {return fedata;}
		std::vector<double> getfmudata() {return fmudata;}
		std::vector<double> getchi2data() {return chi2data;}
		std::string getchi2file() {return schi2file;}	
	private:
		std::vector<double> fedata, fmudata, chi2data;
		std::string schi2file = "";
		std::string year;
		InterpMultilinear<2, double> *interp2d;		
};

class neutrinodecay : public flavorregion
{
public:
	neutrinodecay():flavorregion(){};
	~neutrinodecay(){};
	double rho(double z);
	double h(double z);
	double Z2(double z);
	void calcmatrix(std::vector<double> oscinput);
	double Delta(double kappa, double EnuTeV, double z);
	double D(double kappa, double EnuTeV, double z);
	void cbeta(std::vector<double> comp_i, std::vector<std::vector<double>> &cmatrix);
	double Eint(double c1, double c2, double c3, double Delta2, double Delta3, double EnuTeV, double z);
	double trapz(std::function<double(double)> func, double start, double end, int Nbins);
	void flavorcomp(double kappa2, double kappa3, std::vector<double> comp_i, std::vector<double> oscinput, std::vector<double> &comp_f);
private:
	//parameters for neutrino decay
	const double a = 1.665, b = 1 - a, c = 1.425;
	const double LH = 4.4532; //Gpc
	const double prefactor = 1.02867e5; //for unit conversion
	//parameters for AGN profile
	const double zmax = 4., zc = 1.5, m = 1.5;
	const double Omegam = 0.3158, OmegaL = 1-Omegam; 
	//parameters for neutrino flux
	//const double gamma = 2.87;
	const double gamma = 2.50;
	const double Enumin = 60.; //in TeV
	const double Enumax = 1.e4;
	//for mixing matrix
	std::vector<std::vector<double>> Vsq;


	
};

/* /////////////////////////////////////////////////////////////////////////////

DATA FROM INDIVIDUAL EXPERIMENTS BELOW

 ///////////////////////////////////////////////////////////////////////////// */

/* Data from the JUNO experiment. Numbers taken from xxxx.xxxxx p.x Tab. y  */
class JUNO : public oscillationexperiment
{
	public:
		JUNO();
		~JUNO(){};
		void setosc(oscillationparams &osc);
	private:
		double t12sqsigma = 0.002, dm21sqsigma = 0.044e-5, dm31sqsigma = 0.017e-3;
};




class DUNE : public oscillationexperiment
{
	public:
		DUNE(const std::string ordering = "NO", const std::string t23oct = "upper", const double dcpbest = 1.);
		~DUNE(){};
		void setosc(oscillationparams &osc);
	private:
		double t13sqsigma = 0.0021, dm31sqsigma = 0.0175e-3;

};

/* Data from the HyperK experiment. Numbers from https://arxiv.org/pdf/1805.04163.pdf */
class HYPERK : public oscillationexperiment
{
	public:
		HYPERK(const std::string ordering = "NO", const std::string t23oct = "upper", const double dcpbest = 1.);
		~HYPERK(){};
		void setosc(oscillationparams &osc);
	private:
	    double dcpsigma = 23.; // degrees
	    double t23sqsigma;
	    double dm23sigma = 1.25e-5;
    //double get_ds23(oscillationparams &osc);
};


class NUFIT : public oscillationexperiment
{
	public:
		NUFIT(const std::string ordering = "NO", const double nufitversion = 5.0);
		~NUFIT(){};
		void setosc(oscillationparams &osc);

};

class prior
{
	public:
		prior();
		~prior(){};
		double rand01();
		double flatPrior(double r, double x1, double x2);
		std::vector<double> randomInitialFlavor();
		std::vector<double> random2d();
		//double gaussianPrior(double r, double x1, double x2);
	private:
		std::mt19937 *generator;
		dirichlet_distribution<std::mt19937> *distgenerator_diric;
		dirichlet_distribution<std::mt19937> *distgenerator_diric2d;
		std::uniform_real_distribution<double> *distgenerator;
};

class gaussianprior
{
	public:
		gaussianprior(double mean, double std);
		~gaussianprior(){};
		double gaussian();
	private:
		//std::mt19937 *generator;
		boost::variate_generator<std::mt19937, boost::normal_distribution<double> > *distgenerator;
};