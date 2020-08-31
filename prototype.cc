#include "prototype.h"


oscillationparams::oscillationparams(const std::string ordering) : ordering(ordering)
{
	//initialize oscillation parameters with Nufit 5.0, with SK
	if (ordering == "NO")
	{
// Central value, upper error bar, lower error bar
		set12(0.304, 0.012, 0.012);             // sin^2 \theta_12
		set13(0.02219, 0.00062, 0.00063);       // sin^2 \theta_13
		set23(0.573, 0.016, 0.020);             // sin^2 \theta 23
		setdcp(197, 27, 24);                    //\delta_{CP} in degrees
		setdm21(7.42e-5, 0.21e-5, 0.20e-5);     //\Delta m_{21}^2 in eV^2
		setdm31(2.517e-3, 0.026e-3, 0.028e-3);  //\Delta m_{31}^2 in eV^2
	}
	else if (ordering == "IO")
	{
		set12(0.304, 0.013, 0.012);
		set13(0.022838, 0.00063, 0.00062);
		set23(0.575, 0.016, 0.019);
		setdcp(282, 26, 30);
		setdm21(7.42e-5, 0.21e-5, 0.20e-5);
		setdm31(-2.498e-3, 0.028e-3, 0.028e-3);
	}
	else
	{
		std::cout << "You must specify a correct mass ordering, choose NO, IO or leave blank\n" << std::endl;
		exit(1);
	}
}

void oscillationparams::set12(double t12sqb, double t12sqsigmap, double t12sqsigmam)
{
	t12sqbest = t12sqb;
	if (!std::isnan(t12sqsigmap)) t12sqsigmaplus = t12sqsigmap;
	if (!std::isnan(t12sqsigmam)) t12sqsigmaminus = t12sqsigmam;
}
void oscillationparams::set13(double t13sqb, double t13sqsigmap, double t13sqsigmam)
{
	t13sqbest = t13sqb;
	if (!std::isnan(t13sqsigmap)) t13sqsigmaplus = t13sqsigmap;
	if (!std::isnan(t13sqsigmam)) t13sqsigmaminus = t13sqsigmam;
}
void oscillationparams::set23(double t23sqb, double t23sqsigmap, double t23sqsigmam)
{
	t23sqbest = t23sqb;
	if (!std::isnan(t23sqsigmap)) t23sqsigmaplus = t23sqsigmap;
	if (!std::isnan(t23sqsigmam)) t23sqsigmaminus = t23sqsigmam;
}
void oscillationparams::setdm21(double dm21sqb, double dm21sqsigmap, double dm21sqsigmam)
{
	dm21sqbest = dm21sqb;
	if (!std::isnan(dm21sqsigmap)) dm21sqsigmaplus = dm21sqsigmap;
	if (!std::isnan(dm21sqsigmam)) dm21sqsigmaminus = dm21sqsigmam;
}
void oscillationparams::setdm31(double dm31sqb, double dm31sqsigmap, double dm31sqsigmam)
{
	dm31sqbest = dm31sqb;
	if (!std::isnan(dm31sqsigmap)) dm31sqsigmaplus = dm31sqsigmap;
	if (!std::isnan(dm31sqsigmam)) dm31sqsigmaminus = dm31sqsigmam;
}
void oscillationparams::setdcp(double dcpb, double dcpsigmap, double dcpsigmam)
{
	dcpbest = dcpb;
	if (!std::isnan(dcpsigmap)) dcpsigmaplus = dcpsigmap;
	if (!std::isnan(dcpsigmam)) dcpsigmaminus = dcpsigmam;
}

flavorregion::flavorregion() : comp_f(3, 0.){}
double flavorregion::sq(double x)
{
	return x*x;
}
double flavorregion::V11s(double q12, double q13)
{
	return sq(cos(q12)) * sq(cos(q13));
}
double flavorregion::V12s(double q12, double q13)
{
	return sq(sin(q12)) * sq(cos(q13));
}
double flavorregion::V13s(double q13)
{
	return sq(sin(q13));
}
double flavorregion::V21s(double q12, double q13, double q23, double dcp)
{
	return sq(sin(q12)) * sq(cos(q23)) + sq(cos(q12)) * sq(sin(q13)) * sq(sin(q23)) + 2 * sin(q12) * cos(q12) * sin(q23) * cos(q23) *sin(q13) * cos(dcp);
}
double flavorregion::V22s(double q12, double q13, double q23, double dcp)
{
	return sq(cos(q12)) * sq(cos(q23)) + sq(sin(q12)) * sq(sin(q13)) * sq(sin(q23)) - 2 * sin(q12) * cos(q12) * sin(q23) * cos(q23) *sin(q13) * cos(dcp);
}
double flavorregion::V23s(double q13, double q23)
{
	return sq(cos(q13)) * sq(sin(q23));
}
double flavorregion::V31s(double q12, double q13, double q23, double dcp)
{
	return sq(sin(q12)) * sq(sin(q23)) + sq(cos(q12)) * sq(sin(q13)) * sq(cos(q23)) - 2 * sin(q12) * cos(q12) * sin(q23) * cos(q23) *sin(q13) * cos(dcp);
}
double flavorregion::V32s(double q12, double q13, double q23, double dcp)
{
	return sq(cos(q12)) * sq(sin(q23)) + sq(sin(q12)) * sq(sin(q13)) * sq(cos(q23)) + 2 * sin(q12) * cos(q12) * sin(q23) * cos(q23) *sin(q13) * cos(dcp);
}
double flavorregion::V33s(double q13, double q23)
{
	return sq(cos(q13)) * sq(cos(q23));
}
std::vector<double> flavorregion::evolveflavor(std::vector<double> comp_i, std::vector<double> oscinput)
{
	std:vector<std::vector<double>> Vsq(3, std::vector<double>(3, 0.));
	double q12=asin(sqrt(oscinput[0]));
	double q13=asin(sqrt(oscinput[1]));
	double q23=asin(sqrt(oscinput[2]));
	double dcp=oscinput[3]/180.*M_PI;
	Vsq[0][0] = V11s(q12, q13);
	Vsq[0][1] = V12s(q12, q13);
	Vsq[0][2] = V13s(q13);
	Vsq[1][0] = V21s(q12, q13, q23, dcp);
	Vsq[1][1] = V22s(q12, q13, q23, dcp);
	Vsq[1][2] = V23s(q13, q23);
	Vsq[2][0] = V31s(q12, q13, q23, dcp);
	Vsq[2][1] = V32s(q12, q13, q23, dcp);
	Vsq[2][2] = V33s(q13, q23);

	for(int k = 0; k < 3; k++){
		comp_f[k] = 0.;
		for(int i = 0; i < 3; i++){
			double P = 0;
			for(int j = 0; j < 3; j++)
				P += Vsq[i][j] * Vsq[k][j];
			comp_f[k] += P * comp_i[i];
		}
	}
	return comp_f;
}


oscillationexperiment::oscillationexperiment(const std::string chi2file) : schi2file(chi2file){}
void oscillationexperiment::readchi2(std::string chi2file)
{
	std::ifstream file(chi2file);
	if (!file) {
	    std::cerr << "Unable to open file " << chi2file << ", file not found." << std::endl;
	    exit(1);   // call system to stop
	}			
	std::string line;
	while(getline(file, line))
	{
		//if start with # or empty do nothing
		if (!(line.rfind("#", 0) == 0) && (line.find_first_not_of(" ") != std::string::npos))
		{
			double tmp1, tmp2, tmp3;
			sscanf(line.c_str(), "%lf %lf %lf", &tmp1, &tmp2, &tmp3);
			t23sqdata.push_back(tmp1);
			dcpdata.push_back(tmp2);
			chi2data.push_back(tmp3);
		}
	}
	t23sqdata.erase(remove_duplicates(t23sqdata.begin(), t23sqdata.end()), t23sqdata.end());
	dcpdata.erase(remove_duplicates(dcpdata.begin(), dcpdata.end()), dcpdata.end());
	//normalize chi2
	double minchi2 = *min_element(chi2data.begin(), chi2data.end());
	for(int i = 0; i < chi2data.size(); i++) chi2data[i] -= minchi2;
}


likelihood::likelihood(experimentlist explist, oscillationparams &oscparams)
{
	oscexperiments(explist, oscparams);
	osc = oscparams;
}
void likelihood::oscexperiments(experimentlist explist, oscillationparams &osc)
{
	for(int i = 0; i < explist.size(); i++)
	{
		explist[i]->setosc(osc);
		if (explist[i]->getchi2file() != "")
		{
			t23sqdata = explist[i]->gett23sqdata();
			dcpdata = explist[i]->getdcpdata();
			chi2data = explist[i]->getchi2data();
			schi2file = explist[i]->getchi2file();
		}
	}
	if (schi2file != "")
	{
		std::vector< std::vector<double>::iterator > grid_iter_list;
		grid_iter_list.push_back(t23sqdata.begin());
		grid_iter_list.push_back(dcpdata.begin());
		array<int,2> grid_sizes;
		grid_sizes[0] = t23sqdata.size();
		grid_sizes[1] = dcpdata.size();
		auto interp_ML = new InterpMultilinear<2, double> (grid_iter_list.begin(), grid_sizes.begin(), chi2data.data(), chi2data.data() + chi2data.size());
		interp2d = interp_ML;
	}
}
double likelihood::sq(double x)
{
	return x*x;
}		
double likelihood::chisq(std::vector<double> oscinput)
{

	double t12s=oscinput[0];
	double t13s=oscinput[1];
	double t23s=oscinput[2];
	double dcp=oscinput[3];
    return sq(t12s-osc.get12best())/sq((t12s>osc.get12best())?osc.get12plus():osc.get12minus())
    +sq(t13s-osc.get13best())/sq((t13s>osc.get13best())?osc.get13plus():osc.get13minus())
    +sq(t23s-osc.get23best())/sq((t23s>osc.get23best())?osc.get23plus():osc.get23minus())
    +sq(dcp-osc.getdcpbest())/sq((dcp>osc.getdcpbest())?osc.getdcpplus():osc.getdcpminus());
}

double likelihood::chisq1213(std::vector<double> oscinput)
{
	double t12s=oscinput[0];
	double t13s=oscinput[1];
    return sq(t12s-osc.get12best())/sq((t12s>osc.get12best())?osc.get12plus():osc.get12minus())
    +sq(t13s-osc.get13best())/sq((t13s>osc.get13best())?osc.get13plus():osc.get13minus());
}

double likelihood::chisq23dcp(std::vector<double> oscinput)
{
	double t23s=oscinput[2];
	double dcp=oscinput[3];
	array<double,2> args = {t23s, dcp};
	if (schi2file != "")  return interp2d->interp(args.begin());
	else return 0.;

}
double likelihood::chisqfromdata(std::vector<double> oscinput)
{
	return chisq1213(oscinput) + chisq23dcp(oscinput);
}		


/* /////////////////////////////////////////////////////////////////////////////

DATA FROM INDIVIDUAL EXPERIMENTS BELOW

 ///////////////////////////////////////////////////////////////////////////// */

/* Data from the JUNO experiment. Numbers taken from xxxx.xxxxx p.x Tab. y  */

JUNO::JUNO() : oscillationexperiment() {}
void JUNO::setosc(oscillationparams &osc)
{
	//set the standard deviation to be the smaller one between the curret value and
	//the expected experimental value to be considered
	if (t12sqsigma < osc.get12plus()) osc.set12(osc.get12best(), t12sqsigma, osc.get12minus());
	if (t12sqsigma < osc.get12minus()) osc.set12(osc.get12best(), osc.get12plus(), t12sqsigma);
	if (dm21sqsigma < osc.getdm21plus()) osc.setdm21(osc.getdm21best(), dm21sqsigma, osc.getdm21minus());
	if (dm21sqsigma < osc.getdm21minus()) osc.setdm21(osc.getdm21best(), osc.getdm21plus(), dm21sqsigma);
	if (dm31sqsigma < osc.getdm31plus()) osc.setdm31(osc.getdm31best(), dm31sqsigma, osc.getdm31minus());
	if (dm31sqsigma < osc.getdm31minus()) osc.setdm31(osc.getdm31best(), osc.getdm31plus(), dm31sqsigma);
}




DUNE::DUNE(std::string chi2file) : oscillationexperiment(chi2file)
{
	readchi2(chi2file);
}
void DUNE::setosc(oscillationparams &osc) 
{
	if (t13sqsigma < osc.get13plus()) osc.set13(osc.get13best(), t13sqsigma, osc.get13minus());
	if (t13sqsigma < osc.get13minus()) osc.set13(osc.get13best(), osc.get13plus(), t13sqsigma);
	if (dm31sqsigma < osc.getdm31plus()) osc.setdm31(osc.getdm31best(), dm31sqsigma, osc.getdm31minus());
	if (dm31sqsigma < osc.getdm31minus()) osc.setdm31(osc.getdm31best(), osc.getdm31plus(), dm31sqsigma);			
}


/* Data from the HyperK experiment. Numbers from https://arxiv.org/pdf/1805.04163.pdf */

HYPERK::HYPERK(std::string chi2file) : oscillationexperiment(chi2file)
{
	readchi2(chi2file);
}
void HYPERK::setosc(oscillationparams &osc) {}
//determination of deltaCP from hyperK is still premature
/*
void HYPERK::setosc(oscillationparams &osc)
{
	//set the standard deviation to be the smaller one between the curret value and
	//the expected experimental value to be considered
	if (dcpsigma < osc.getdcpplus()) osc.setdcp(osc.getdcpbest(), dcpsigma, osc.getdcpminus());
    if (dcpsigma < osc.getdcpminus()) osc.setdcp(osc.getdcpbest(), osc.getdcpplus(),dcpsigma);
    t23sqsigma = get_ds23(osc);
    if (t23sqsigma < osc.get23plus()) osc.set23(osc.get23best(), t23sqsigma, osc.get23minus());
    if (t23sqsigma < osc.get23minus()) osc.set23(osc.get23best(), osc.get23plus(), t23sqsigma);

    //Delta m23
    if (dm23sigma < osc.getdm31minus()) osc.setdm31(osc.getdm31best(), dm23sigma, osc.getdm31minus());
    if (dm23sigma < osc.get23minus()) osc.set23(osc.get23best(), osc.getdm31plus(), dm23sigma);

}

double HYPERK::get_ds23(oscillationparams &osc)
{
	double dst;
	double sin23 = osc.get12best();
	// Polynomial fit to table in TDR. This is wrong
	dst =  -3.8*sin23*sin23 + 3.83*sin23 - 0.948;
	if (dst < 0.0)
	{
	std::cout << "This fit for hyperK was a great idea but it doesn't work, sorry";
	dst  = 0.0;
	}
    return dst;
}
*/

NUFIT::NUFIT(std::string chi2file) : oscillationexperiment(chi2file)
{
	readchi2(chi2file);
}
void NUFIT::setosc(oscillationparams &osc) {}


double prior::rand01()
{
	return rand() / double(RAND_MAX);
}
double prior::flatPrior(double r, double x1, double x2)
{
	return x1 + r * ( x2 - x1 );
}
/*
double prior::gaussianPrior(double r, double x1, double x2)
{
	if( x2 <= 0.0 )
	{
		printf("sigma <= 0 in routine gaussianPrior\n");
		exit(1);
	}
	
	if( r <= 1.E-16 || ( 1.0 - r ) <= 1.E-16 )
       	return ( -1.0E32 );
    else
       	return ( x1 + x2 * sqrt( 2.0 ) * dierfc( 2.0 * ( 1.0 - r ) ) );
}
*/