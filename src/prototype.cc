#include "prototype.h"

//print an m * n matrix from a 1d vector
void printmatrix(std::vector<double> v, int m, int n)
{
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			printf("%1.6lf ", v[i * m + j]);
		}
		printf("\n");
	}
}

oscillationparams::oscillationparams(const std::string ordering, const double NUFITversion) : ordering(ordering), version(NUFITversion)
{
	if (ordering == "NO")
	{
		// Central value, upper error bar, lower error bar
		switch(int(version*10))
		{
			//initialize oscillation parameters with Nufit 5.0, with SK
			case 50:
				set12(0.304, 0.012, 0.012);             // sin^2 \theta_12
				set13(0.02219, 0.00062, 0.00063);       // sin^2 \theta_13
				set23(0.573, 0.016, 0.020);             // sin^2 \theta 23
				setdcp(197., 27., 24.);                    //\delta_{CP} in degrees
				setdm21(7.42e-5, 0.21e-5, 0.20e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.517e-3, 0.026e-3, 0.028e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 4.1, with SK
			case 41:
				set12(0.310, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.02237, 0.00066, 0.00065);       // sin^2 \theta_13
				set23(0.563, 0.018, 0.024);             // sin^2 \theta 23
				setdcp(221., 39., 28.);                    //\delta_{CP} in degrees
				setdm21(7.39e-5, 0.21e-5, 0.20e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.528e-3, 0.029e-3, 0.031e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 4.0, with SK
			case 40:
				set12(0.310, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.02240, 0.00065, 0.00066);       // sin^2 \theta_13
				set23(0.582, 0.015, 0.019);             // sin^2 \theta 23
				setdcp(217., 40., 28.);                    //\delta_{CP} in degrees
				setdm21(7.39e-5, 0.21e-5, 0.20e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.525e-3, 0.033e-3, 0.031e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 3.2
			case 32:
				set12(0.307, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.02206, 0.00075, 0.00075);       // sin^2 \theta_13
				set23(0.538, 0.033, 0.069);             // sin^2 \theta 23
				setdcp(234., 43., 31.);                    //\delta_{CP} in degrees
				setdm21(7.40e-5, 0.21e-5, 0.20e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.494e-3, 0.033e-3, 0.031e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 3.1
			case 31:
				set12(0.307, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.02195, 0.00075, 0.00074);       // sin^2 \theta_13
				set23(0.565, 0.025, 0.120);             // sin^2 \theta 23
				setdcp(228., 51., 33.);                    //\delta_{CP} in degrees
				setdm21(7.40e-5, 0.21e-5, 0.20e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.515e-3, 0.035e-3, 0.035e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 3.0
			case 30:
				set12(0.306, 0.012, 0.012);             // sin^2 \theta_12
				set13(0.02166, 0.00075, 0.00075);       // sin^2 \theta_13
				set23(0.441, 0.027, 0.021);             // sin^2 \theta 23
				setdcp(261., 51., 59.);                    //\delta_{CP} in degrees
				setdm21(7.50e-5, 0.19e-5, 0.17e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.524e-3, 0.039e-3, 0.040e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with 2.2
			case 22:
				set12(0.308, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.02163, 0.00074, 0.00074);       // sin^2 \theta_13
				set23(0.440, 0.023, 0.019);             // sin^2 \theta 23
				setdcp(289., 38., 51);                    //\delta_{CP} in degrees
				setdm21(7.49e-5, 0.19e-5, 0.17e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.526e-3, 0.039e-3, 0.037e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 2.1 "LID"
			case 21:
				set12(0.308, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.0219, 0.0010, 0.0010);       // sin^2 \theta_13
				set23(0.451, 0.038, 0.025);             // sin^2 \theta 23
				setdcp(303., 39., 50.);                    //\delta_{CP} in degrees
				setdm21(7.49e-5, 0.19e-5, 0.17e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.477e-3, 0.042e-3, 0.042e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 2.0
			case 20:
				set12(0.304, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.0218, 0.0010, 0.0010);       // sin^2 \theta_13
				set23(0.452, 0.052, 0.028);             // sin^2 \theta 23
				setdcp(306., 39., 70.);                    //\delta_{CP} in degrees
				setdm21(7.50e-5, 0.19e-5, 0.17e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.457e-3, 0.047e-3, 0.047e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 1.3, lower octant with RSBL
			case 13:
				set12(0.304, 0.012, 0.012);             // sin^2 \theta_12
				set13(0.0219, 0.0010, 0.0011);       // sin^2 \theta_13
				set23(0.451, 0.001, 0.001);             // sin^2 \theta 23
				setdcp(251., 67., 59.);                    //\delta_{CP} in degrees
				setdm21(7.50e-5, 0.19e-5, 0.17e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.458e-3, 0.002e-3, 0.002e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 1.2, lower octant with RSBL
			case 12:
				set12(0.306, 0.012, 0.012);             // sin^2 \theta_12
				set13(0.0231, 0.0019, 0.0019);       // sin^2 \theta_13
				set23(0.446, 0.008, 0.008);             // sin^2 \theta 23
				setdcp(266., 55., 63.);                    //\delta_{CP} in degrees
				setdm21(7.45e-5, 0.19e-5, 0.16e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.417e-3, 0.014e-3, 0.014e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 1.1 with RSBL
			case 11:
				set12(0.306, 0.012, 0.012);             // sin^2 \theta_12
				set13(0.0231, 0.0023, 0.022);       // sin^2 \theta_13
				set23(0.437, 0.061, 0.031);             // sin^2 \theta 23
				setdcp(341., 58., 46.);                    //\delta_{CP} in degrees
				setdm21(7.45e-5, 0.19e-5, 0.16e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.421e-3, 0.022e-3, 0.023e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 1.0, lower octant with RSBL
			case 10:
				set12(0.302, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.0227, 0.0023, 0.0024);       // sin^2 \theta_13
				set23(0.413, 0.037, 0.025);             // sin^2 \theta 23
				setdcp(300., 66., 138.);                    //\delta_{CP} in degrees
				setdm21(7.50e-5, 0.18e-5, 0.19e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.473e-3, 0.070e-3, 0.067e-3);  //\Delta m_{31}^2 in eV^2
				break;
			default:
				std::cout << "Nufit version doesn't match any record, exit." << std::endl;
				exit(1);																																																									
		}

	}
	else if (ordering == "IO")
	{
		// Central value, upper error bar, lower error bar
		switch(int(version*10))
		{
			//initialize oscillation parameters with Nufit 5.0, with SK
			case 50:
				set12(0.304, 0.013, 0.012);
				set13(0.022838, 0.00063, 0.00062);
				set23(0.575, 0.016, 0.019);
				setdcp(282, 26, 30);
				setdm21(7.42e-5, 0.21e-5, 0.20e-5);
				setdm31(-2.498e-3, 0.028e-3, 0.028e-3);
				break;
			//initialize oscillation parameters with Nufit 4.1, with SK
			case 41:
				set12(0.310, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.02259, 0.00065, 0.00065);       // sin^2 \theta_13
				set23(0.565, 0.017, 0.022);             // sin^2 \theta 23
				setdcp(282., 23., 25.);                    //\delta_{CP} in degrees
				setdm21(7.39e-5, 0.21e-5, 0.20e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.510e-3, 0.030e-3, 0.031e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 4.0, with SK
			case 40:
				set12(0.310, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.02263, 0.00065, 0.00066);       // sin^2 \theta_13
				set23(0.582, 0.015, 0.018);             // sin^2 \theta 23
				setdcp(280., 25., 28.);                    //\delta_{CP} in degrees
				setdm21(7.39e-5, 0.21e-5, 0.20e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(2.515e-3, 0.034e-3, 0.031e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 3.2
			case 32:
				set12(0.307, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.02227, 0.00074, 0.00074);       // sin^2 \theta_13
				set23(0.554, 0.023, 0.033);             // sin^2 \theta 23
				setdcp(278., 26., 29.);                    //\delta_{CP} in degrees
				setdm21(7.40e-5, 0.21e-5, 0.20e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(-2.465e-3, 0.032e-3, 0.031e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 3.1
			case 31:
				set12(0.307, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.02212, 0.00074, 0.00073);       // sin^2 \theta_13
				set23(0.572, 0.021, 0.028);             // sin^2 \theta 23
				setdcp(281., 30., 33.);                    //\delta_{CP} in degrees
				setdm21(7.40e-5, 0.21e-5, 0.20e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(-2.483e-3, 0.034e-3, 0.035e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 3.0
			case 30:
				set12(0.306, 0.012, 0.012);             // sin^2 \theta_12
				set13(0.02179, 0.00076, 0.00076);       // sin^2 \theta_13
				set23(0.587, 0.020, 0.024);             // sin^2 \theta 23
				setdcp(277., 40., 46.);                    //\delta_{CP} in degrees
				setdm21(7.50e-5, 0.19e-5, 0.17e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(-2.514e-3, 0.038e-3, 0.041e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with 2.2
			case 22:
				set12(0.308, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.02175, 0.00075, 0.00074);       // sin^2 \theta_13
				set23(0.584, 0.018, 0.022);             // sin^2 \theta 23
				setdcp(269., 39., 45.);                    //\delta_{CP} in degrees
				setdm21(7.49e-5, 0.19e-5, 0.17e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(-2.518e-3, 0.038e-3, 0.037e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 2.1 "LID"
			case 21:
				set12(0.308, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.0219, 0.0010, 0.0010);       // sin^2 \theta_13
				set23(0.576, 0.023, 0.033);             // sin^2 \theta 23
				setdcp(262., 51., 57.);                    //\delta_{CP} in degrees
				setdm21(7.49e-5, 0.19e-5, 0.17e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(-2.465e-3, 0.041e-3, 0.043e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 2.0
			case 20:
				set12(0.304, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.0219, 0.0011, 0.0010);       // sin^2 \theta_13
				set23(0.579, 0.025, 0.037);             // sin^2 \theta 23
				setdcp(254., 63., 62.);                    //\delta_{CP} in degrees
				setdm21(7.50e-5, 0.19e-5, 0.17e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(-2.449e-3, 0.048e-3, 0.047e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 1.3, upper octant with RSBL
			case 13:
				set12(0.304, 0.012, 0.012);             // sin^2 \theta_12
				set13(0.0219, 0.0010, 0.0011);       // sin^2 \theta_13
				set23(0.577, 0.027, 0.0035);             // sin^2 \theta 23
				setdcp(251., 67., 59.);                    //\delta_{CP} in degrees
				setdm21(7.50e-5, 0.19e-5, 0.17e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(-2.448e-3, 0.047e-3, 0.047e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 1.2, upper octant with RSBL
			case 12:
				set12(0.306, 0.012, 0.012);             // sin^2 \theta_12
				set13(0.0231, 0.0019, 0.0019);       // sin^2 \theta_13
				set23(0.593, 0.027, 0.043);             // sin^2 \theta 23
				setdcp(266., 55., 63.);                    //\delta_{CP} in degrees
				setdm21(7.45e-5, 0.19e-5, 0.16e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(-2.411e-3, 0.062e-3, 0.062e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 1.1 with RSBL
			case 11:
				set12(0.306, 0.012, 0.012);             // sin^2 \theta_12
				set13(0.0231, 0.0023, 0.022);       // sin^2 \theta_13
				set23(0.437, 0.061, 0.031);             // sin^2 \theta 23
				setdcp(341., 58., 46.);                    //\delta_{CP} in degrees
				setdm21(7.45e-5, 0.19e-5, 0.16e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(-2.410e-3, 0.062e-3, 0.063e-3);  //\Delta m_{31}^2 in eV^2
				break;
			//initialize oscillation parameters with Nufit 1.0, upper octant with RSBL
			case 10:
				set12(0.302, 0.013, 0.012);             // sin^2 \theta_12
				set13(0.0227, 0.0023, 0.0024);       // sin^2 \theta_13
				set23(0.594, 0.021, 0.022);             // sin^2 \theta 23
				setdcp(300., 66., 138.);                    //\delta_{CP} in degrees
				setdm21(7.50e-5, 0.18e-5, 0.19e-5);     //\Delta m_{21}^2 in eV^2
				setdm31(-2.427e-3, 0.042e-3, 0.065e-3);  //\Delta m_{31}^2 in eV^2
				break;
			default:
				std::cout << "Nufit version doesn't match any record, exit." << std::endl;
				exit(1);																																																									
		}		
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


nonunitflavorregion::nonunitflavorregion(const int year, const std::string oscoption)
{
	fillchi2(year, oscoption);
}

void nonunitflavorregion::readchi2(std::string chi2file)
{
	std::ifstream file(chi2file);
	if (!file) {
	    std::cerr << "Unable to open file " << chi2file << ", file not found." << std::endl;
	    exit(1);   // call system to stop
	}			
	std::string line;
	matrixdata U;
	while(getline(file, line))
	{
		//if start with # or empty do nothing
		if (!(line.rfind("#", 0) == 0) && (line.find_first_not_of(" ") != std::string::npos))
		{
			double tmp1, tmp2;
			sscanf(line.c_str(), "%lf %lf", &tmp1, &tmp2);
			U.Usqdata.push_back(tmp1);
			U.Usqchi2.push_back(tmp2);
		}
	}
	//normalize chi2
	double minchi2 = *min_element(U.Usqchi2.begin(), U.Usqchi2.end());
	for(int i = 0; i < U.Usqchi2.size(); i++) U.Usqchi2[i] -= minchi2;
	//cut too large chi2
	auto it = U.Usqchi2.begin();
	auto itd = U.Usqdata.begin();
	while (it != U.Usqchi2.end())
	{
		if (*it > 16.)
		{
			itd = U.Usqdata.erase(itd);
			it = U.Usqchi2.erase(it);
		}
		else {++it; ++itd;}
	}
	Usq.push_back(U);
}

void nonunitflavorregion::fillchi2(const int year, const std::string oscoption)
{
	std::string syear;
	if (year == 2020) syear = "Current";
	else if (year == 2040) syear = "Future";
	else std::cerr << "Wrong year specified" << ", choose 2020 or 2040." << std::endl;
	char flavor[] = "emt";
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			std::string fname = "data/non-unitarity_Chi2/U";
			if (oscoption == "agnostic") fname += flavor[i]+std::to_string(j+1)+"Sq_"+syear+"_Chi2.dat";
			else if (oscoption == "submatrix") fname += flavor[i]+std::to_string(j+1)+"Sq_"+syear+"_Submatrix"+"_Chi2.dat";
			else std::cerr << "Wrong assumption specified" << ", choose submatrix or agnostic." << std::endl;
			readchi2(fname);
			std::vector< std::vector<double>::iterator > grid_iter_list;
			grid_iter_list.push_back(Usq[i*3+j].Usqdata.begin());
			array<int,1> grid_sizes;
			grid_sizes[0] = Usq[i*3+j].Usqdata.size();
			auto interp_ML = new InterpMultilinear<1, double> (grid_iter_list.begin(), grid_sizes.begin(), Usq[i*3+j].Usqchi2.data(), Usq[i*3+j].Usqchi2.data() + Usq[i*3+j].Usqchi2.size());
			interp1d[i*3+j] = interp_ML;

		}
	}

	//find the bestfit of each matrix element
	for (int i = 0; i < Usq.size(); ++i)
	{
		int ind = std::min_element(Usq[i].Usqchi2.begin(),Usq[i].Usqchi2.end()) - Usq[i].Usqchi2.begin();
		Usqbest.push_back(Usq[i].Usqdata[ind]);
	}
}

double nonunitflavorregion::chisq(std::vector<double> Usqinput)
{
	double chi2 = 0.;
	for (int i = 0; i < 9; ++i)
	{
		array<double,1> args = {Usqinput[i]};
		chi2 += interp1d[i]->interp(args.begin());
	}
	return chi2;
}

std::vector<double> nonunitflavorregion::evolvefromflavor(std::vector<double> comp_i, std::vector<double> Usqinput, bool normalized)
{
	std:vector<std::vector<double>> Vsq(3, std::vector<double>(3, 0.));
	std::vector<double> Nalpha (3, 0.);
	for (int i = 0; i < 9; ++i)
	{
		Vsq[i/3][i%3] = Usqinput[i];
	}

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			Nalpha[i] += Vsq[i][j];
		}
	}

	for(int k = 0; k < 3; k++){
		comp_f[k] = 0.;
		for(int i = 0; i < 3; i++){
			double P = 0;
			for(int j = 0; j < 3; j++)
				P += Vsq[i][j] * Vsq[k][j] / Nalpha[i] / Nalpha[k];
			comp_f[k] += P * comp_i[i];
		}
	}
	if (normalized == true)
	{
		//sum of final flavor composition
		double fsum = accumulate(comp_f.begin(), comp_f.end(), 0.);
		for (int i = 0; i < comp_f.size(); ++i) comp_f[i] /= fsum;
	}
	return comp_f;	
}


//flavorregion::flavorregion() : comp_f(3, 0.){}
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
std::vector<double> flavorregion::evolvefromflavor(std::vector<double> comp_i, std::vector<double> oscinput)
{
	std::vector<double> comp_f {std::vector<double>(3,0.)};
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

std::vector<double> flavorregion::evolvebackfromflavor(std::vector<double> comp_f, std::vector<double> oscinput)
{
	std::vector<double> comp_i {std::vector<double>(3,0.)};
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
	//calculate probability matrix
	double P[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			for (int k = 0; k < 3; ++k)
				P[i*3+j] += Vsq[j][k] * Vsq[i][k];
	//calculate P^-1
	double *Pinv = symmetricmatrixinversion(P, 3);
	for	(int i = 0; i < 3; ++i) 
		for (int j = 0; j < 3; ++j)
			comp_i[i] += Pinv[i*3+j] * comp_f[j];
	return comp_i;
}


std::vector<double> flavorregion::evolvefrommass(std::vector<double> comp_i, std::vector<double> oscinput)
{
	std::vector<double> comp_f {std::vector<double>(3,0.)};
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
			comp_f[k] += Vsq[k][i] * comp_i[i];
		}
	}
	return comp_f;
}

//Cholesky decomposition, return the lower triangle matrix
double* flavorregion::cholesky(double *A, int n) {
    double *L = (double*)calloc(n * n, sizeof(double));
    if (L == NULL)
        exit(EXIT_FAILURE);
 
    for (int i = 0; i < n; i++)
        for (int j = 0; j < (i+1); j++) {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += L[i * n + k] * L[j * n + k];
            L[i * n + j] = (i == j) ?
                           sqrt(A[i * n + i] - s) :
                           (1.0 / L[j * n + j] * (A[i * n + j] - s));
        }
 
    return L;
}
 

//inversion of lower triangular matrix and return the full inversed matrix
double* flavorregion::lowermatrixinversion(double *A, int n) {
    double *L = (double*)calloc(n * n, sizeof(double));
    double *Lout = (double*)calloc(n * n, sizeof(double));
    if (L == NULL)
        exit(EXIT_FAILURE);
    for (int i = 0; i < n; i++) {
        L[i * n + i] = 1. / A[i * n + i];
        for (int j = 0; j < i; j++) {
            double s = 0.;
            for (int k = j; k < i; k++)
                s += A[i * n + k] * L[k * n + j];
            L[i * n + j] = - s * L[i * n + i];
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                Lout[i * n + j] += L[k * n + i] * L[k * n + j];
    return Lout;
}


//inversion of a symmetric matrix
double* flavorregion::symmetricmatrixinversion(double *A, int n) {
	double *L = cholesky(A, n);
	double *Lout = lowermatrixinversion(L, n);
	return Lout;
}

void flavorregion::show_matrix(double *A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%2.5f ", A[i * n + j]);
        printf("\n");
    }
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
			//if not the first chi2 map
			if ((t23sqdatamin != 0.) || (t23sqdatamax != 0.) || (dcpdatamin != 0.) || (dcpdatamax != 0.))
			{
				std::vector<double> v = explist[i]->gett23sqdata();
				double t23min = *std::min_element(v.begin(), v.end());
				double t23max = *std::max_element(v.begin(), v.end());
				v = explist[i]->getdcpdata();
				double dcpmin = *std::min_element(v.begin(), v.end());
				double dcpmax = *std::max_element(v.begin(), v.end());
				//check if the chi2 maps match, these conditions can be relaxed later
				if (((t23sqdatamin != t23min) || (t23sqdatamax != t23max) || (dcpdatamin != dcpmin) || (dcpdatamax != dcpmax)) ||
					((t23sqdata.size() != (explist[i]->gett23sqdata()).size()) || (dcpdata.size() != (explist[i]->getdcpdata()).size())))
				{
					std::cout << "The data in chi2 maps provided don't match, exit." << std::endl;
					exit(1);
				}
				else
				{
					//sum up chi2 from different maps
					std::transform (chi2data.begin(), chi2data.end(), (explist[i]->getchi2data()).begin(), chi2data.begin(), std::plus<int>());
					schi2file += ", " + explist[i]->getchi2file();
				}

			}
			//if first chi2 map, fill up everything
			else
			{
				t23sqdata = explist[i]->gett23sqdata();
				dcpdata = explist[i]->getdcpdata();
				chi2data = explist[i]->getchi2data();
				schi2file = explist[i]->getchi2file();				
				//t23sqdatamin = *std::min_element((explist[i]->gett23sqdata()).begin(), (explist[i]->gett23sqdata()).end());
				t23sqdatamin = *std::min_element(t23sqdata.begin(), t23sqdata.end());
				t23sqdatamax = *std::max_element(t23sqdata.begin(), t23sqdata.end());
				dcpdatamin = *std::min_element(dcpdata.begin(), dcpdata.end());
				dcpdatamax = *std::max_element(dcpdata.begin(), dcpdata.end());
			}

		}
	}
	if (schi2file != "")
	{
		//normalize chi2 again in case the minima of different chi2 maps don't match
		double minchi2 = *min_element(chi2data.begin(), chi2data.end());
		for(int i = 0; i < chi2data.size(); i++) chi2data[i] -= minchi2;
		//initialize the interpolatoin function		
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


likelihood_ice::likelihood_ice(std::string year) : year(year)
{
	
	if (year != "2015")
	{
		std::string fname = "data/IceCube_Chi2/icecube";
		if (year == "2020") fname += "_8yr.txt";
		else if (year == "2028") fname += "_15yr.txt";
		else if (year == "2040") fname += "+gen2-inice.txt";
		else if (year == "2040-comb") fname = "data/IceCube_Chi2/combineexp.txt";
		else  std::cerr << "Wrong year specified, choose in between 2015, 2020, 2028, 2040 and 2040-comb." << std::endl; 
		readchi2(fname);
		//std::cout << fname << std::endl;

		//initialize the interpolatoin function		
		std::vector< std::vector<double>::iterator > grid_iter_list;
		grid_iter_list.push_back(fedata.begin());
		grid_iter_list.push_back(fmudata.begin());
		array<int,2> grid_sizes;
		grid_sizes[0] = fedata.size();
		grid_sizes[1] = fmudata.size();
		auto interp_ML = new InterpMultilinear<2, double> (grid_iter_list.begin(), grid_sizes.begin(), chi2data.data(), chi2data.data() + chi2data.size());
		interp2d = interp_ML;
	}	

}

void likelihood_ice::readchi2(std::string chi2file)
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
			fedata.push_back(tmp1);
			fmudata.push_back(tmp2);
			chi2data.push_back(tmp3);
		}
	}
	fedata.erase(remove_duplicates(fedata.begin(), fedata.end()), fedata.end());
	fmudata.erase(remove_duplicates(fmudata.begin(), fmudata.end()), fmudata.end());
	//normalize chi2
	double minchi2 = *min_element(chi2data.begin(), chi2data.end());
	for(int i = 0; i < chi2data.size(); i++) chi2data[i] -= minchi2;
}

double likelihood_ice::chisqice(std::vector<double> flavinput)
{
	if (year == "2015") return chisqice_2015(flavinput);
	else
	{
		double fe = flavinput[0];
		double fmu = flavinput[1];
		array<double,2> args = {fe, fmu};
		return interp2d->interp(args.begin());		
	}	
}

double likelihood_ice::chisqice_2015(std::vector<double> flavinput)
{
	double f_eE = flavinput[0];
	double f_mE = flavinput[1];

	double f_eE_bf = 0.50;
	double f_mE_bf = 0.50;
	double sigma_f_eE = 0.970;
	double sigma_f_mE = 0.365;
	double correlation = -0.60;
	double norm = 39.4056;
    double f1 = 1./(2.0*M_PI*sigma_f_eE*sigma_f_mE*sqrt(1.0-correlation*correlation));
    double f2 = 0.5/(1.0-correlation*correlation);
    double f3 = pow(f_eE-f_eE_bf, 2.0) / (sigma_f_eE*sigma_f_eE);
    double f4 = pow(f_mE-f_mE_bf, 2.0) / (sigma_f_mE*sigma_f_mE);
    double f5 = 2.0*correlation*(f_eE-f_eE_bf)*(f_mE-f_mE_bf)/sigma_f_eE/sigma_f_mE;

    // offset is the 2D Gaussian evaluated at the best-fit point
    double offset = f1;

    return norm*(offset-f1*exp(-f2*(f3+f4+f5)));

}

void neutrinodecay::calcmatrix(std::vector<double> oscinput)
{
	std:vector<std::vector<double>> Vsqtmp(3, std::vector<double>(3, 0.));
	double q12=asin(sqrt(oscinput[0]));
	double q13=asin(sqrt(oscinput[1]));
	double q23=asin(sqrt(oscinput[2]));
	double dcp=oscinput[3]/180.*M_PI;

	Vsqtmp[0][0] = V11s(q12, q13);
	Vsqtmp[0][1] = V12s(q12, q13);
	Vsqtmp[0][2] = V13s(q13);
	Vsqtmp[1][0] = V21s(q12, q13, q23, dcp);
	Vsqtmp[1][1] = V22s(q12, q13, q23, dcp);
	Vsqtmp[1][2] = V23s(q13, q23);
	Vsqtmp[2][0] = V31s(q12, q13, q23, dcp);
	Vsqtmp[2][1] = V32s(q12, q13, q23, dcp);
	Vsqtmp[2][2] = V33s(q13, q23);	
	Vsq = Vsqtmp;

	/*
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			std::cout << Vsq[i][j] << " ";
		}
		std::cout << std::endl;
	}
	*/	
}


double neutrinodecay::rho(double z)
{
	if ((z >= 0.) && (z < zc)) return pow((1.+z),m);
	else if ((z >= zc) && (z <= zmax)) return pow((1.+zc),m);
	else {std::cerr << "Error, redshift must be between 0 and " << zmax << std::endl; exit(1);}
}

double neutrinodecay::h(double z)
{
	return sqrt(OmegaL + pow((1.+z), 3)*Omegam);
}

double neutrinodecay::Z2(double z)
{
	return a+b*exp(-c*z);
}

//Delta = kappa*LH*log(Z2)/Enu
double neutrinodecay::Delta(double kappa, double EnuTeV, double z)
{
	return prefactor*kappa*LH/EnuTeV*Z2(z);
}

double neutrinodecay::D(double kappa, double EnuTeV, double z)
{
	return pow(Z2(z), -prefactor*kappa*LH/EnuTeV);
}

//compute c_i^beta
void neutrinodecay::cbeta(std::vector<double> comp_i, std::vector<std::vector<double>> &cmatrix)
{
	double sum = 0.;
	for (int k = 0; k < 3; ++k)
		for (int i = 0; i < 3; ++i)
		{
			cmatrix[i][k] = 0.;
			for (int j = 0; j < 3; ++j)
			{
				cmatrix[i][k] += Vsq[j][i]*comp_i[j];
			}
			cmatrix[i][k] *= Vsq[k][i];
		}
}

double neutrinodecay::Eint(double c1, double c2, double c3, double Delta2, double Delta3, double EnuTeV, double z)
{
	return pow(EnuTeV, 1.-gamma)*(c1/(1.-gamma) + 
		c2*boost::math::tgamma(gamma-1., Delta2)*pow(Delta2, 1.-gamma) + 
		c3*boost::math::tgamma(gamma-1., Delta3)*pow(Delta3, 1.-gamma));
}

double neutrinodecay::trapz(std::function<double(double)> func, double start, double end, int Nbins)
{
    double step = (end-start)/Nbins;
    double integral = 0.5*(func(start)+func(end));
    for (int i = 0; i < Nbins; ++i)
    {
        integral += func(start+step*i);
    }
    integral *= step;
    return integral;
}

void neutrinodecay::flavorcomp(double kappa2, double kappa3, std::vector<double> comp_i, std::vector<double> oscinput, std::vector<double> &comp_f)
{
	calcmatrix(oscinput);
	std:vector<std::vector<double>> cmatrix(3, std::vector<double>(3, 0.));
	cbeta(comp_i, cmatrix);
	double sum = 0;
	for (int i = 0; i < 3; ++i)
	{
		
		auto f = [this, cmatrix, i, kappa2, kappa3](double z) {
			return rho(z)/h(z)*pow((1.+z), -gamma)*
			(Eint(cmatrix[0][i], cmatrix[1][i], cmatrix[2][i], Delta(kappa2, Enumax, z), Delta(kappa3, Enumax, z), Enumax, z)
			-Eint(cmatrix[0][i], cmatrix[1][i], cmatrix[2][i], Delta(kappa2, Enumin, z), Delta(kappa3, Enumin, z), Enumin, z));
		};
		//comp_f[i] = boost::math::quadrature::trapezoidal(f, 0., zmax);
		comp_f[i] = trapz(f, 0., zmax, 100);
		sum += comp_f[i];
	}

	for (int i = 0; i < 3; ++i)
	{
		comp_f[i] /= sum;
	}


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


DUNE::DUNE(const std::string ordering, const std::string t23oct, const double dcpbest)
{
	std::string fname = "data/deltacp_theta23sq_Chi2/deltacp_theta23";
	if (t23oct == "upper") fname += "_";
	else if (t23oct == "lower") fname += "_woct_";
	else if (t23oct == "max") fname += "_max_";
	else {std::cout << "Wrong octant." << std::endl; exit(1);}		
	if (ordering == "NO") fname += "NH";
	else if (ordering == "IO") fname += "IH";
	else {std::cout << "Wrong mass ordering." << std::endl; exit(1);}	
	fname += "_cursol";

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
			
	fname += ".dat";
	schi2file = fname;
	//std::cout << fname << std::endl;
	readchi2(schi2file);	
}

void DUNE::setosc(oscillationparams &osc) 
{
	if (t13sqsigma < osc.get13plus()) osc.set13(osc.get13best(), t13sqsigma, osc.get13minus());
	if (t13sqsigma < osc.get13minus()) osc.set13(osc.get13best(), osc.get13plus(), t13sqsigma);
	if (dm31sqsigma < osc.getdm31plus()) osc.setdm31(osc.getdm31best(), dm31sqsigma, osc.getdm31minus());
	if (dm31sqsigma < osc.getdm31minus()) osc.setdm31(osc.getdm31best(), osc.getdm31plus(), dm31sqsigma);			
}


/* Data from the HyperK experiment. Numbers from https://arxiv.org/pdf/1805.04163.pdf */

HYPERK::HYPERK(const std::string ordering, const std::string t23oct, const double dcpbest)
{
	std:: string fname = "data/deltacp_theta23sq_Chi2/deltacp_theta23";
	if (t23oct == "upper") fname += "_";
	else if (t23oct == "lower") fname += "_woct_";
	else if (t23oct == "max") fname += "_max_";
	else {std::cout << "Wrong octant." << std::endl; exit(1);}		
	if (ordering == "NO") fname += "NH";
	else if (ordering == "IO") fname += "IH";
	else {std::cout << "Wrong mass ordering." << std::endl; exit(1);}
	fname += "_HK";
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
	fname += ".dat";	
	schi2file = fname;
	//std::cout << fname << std::endl;
	readchi2(schi2file);	
}
void HYPERK::setosc(oscillationparams &osc) {}


NUFIT::NUFIT(const std::string ordering, const double nufitversion)
{
	std::string fname = "";
	if (nufitversion >= 2.0)
	{
		fname += "data/deltacp_theta23sq_Chi2/v"+std::to_string(int(nufitversion*10))+".release-data-"+ordering+"_s23sq_dcp.dat";
	}
	if (fname == "") std::cout << "Warning: Nufit chi2 is supposed to be used but no chi2 table served, this could happen is Nufit version < 2.0\n" << std::endl;
	else readchi2(fname);
	schi2file = fname;
	//std::cout << fname << std::endl;
}
void NUFIT::setosc(oscillationparams &osc) {}


prior::prior()
{
	//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	//auto rng = new std::mt19937 (seed);
	auto dis = new std::uniform_real_distribution<double> (0.0, 1.0);
    // initialize rng
    std::random_device rd;
    auto rng = new std::mt19937 (rd());	
	auto dis_diric = new dirichlet_distribution<std::mt19937> ({1, 1, 1});
	auto dis_diric2d = new dirichlet_distribution<std::mt19937> ({1, 1});
	generator = rng;
	distgenerator = dis;
	distgenerator_diric = dis_diric;
	distgenerator_diric2d = dis_diric2d;
}

std::vector<double> prior::randomInitialFlavor()
{
	return (*distgenerator_diric) (*generator);
}

std::vector<double> prior::random2d()
{
	return (*distgenerator_diric2d) (*generator);
}

double prior::rand01()
{
	return (*distgenerator) (*generator);
	//return rand() / double(RAND_MAX);
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


gaussianprior::gaussianprior(double mean, double std)
{
    std::random_device rd;
    auto rng = new std::mt19937 (rd());
    auto gaussiandist = new boost::normal_distribution<double> (mean, std);	
    auto gaussiangen = new boost::variate_generator<std::mt19937, boost::normal_distribution<double> > ((*rng), (*gaussiandist));
	//generator = rng;
	distgenerator = gaussiangen;
}

double gaussianprior::gaussian()
{
	return (*distgenerator) ();
}