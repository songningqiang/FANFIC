#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

class oscillationparams
{
	public:
		oscillationparams(const std::string ordering = "NO") : ordering(ordering)
		{
			//initialize oscillation parameters with Nufit 5.0, with SK
			if (ordering == "NO")
			{
        // Central value, lower error bar, upper error bar
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
		~oscillationparams(){};

		void set12(double t12sqb, double t12sqsigmap = NAN, double t12sqsigmam = NAN)
		{
			t12sqbest = t12sqb;
			if (!std::isnan(t12sqsigmap)) t12sqsigmaplus = t12sqsigmap;
			if (!std::isnan(t12sqsigmam)) t12sqsigmaminus = t12sqsigmam;
		}
		void set13(double t13sqb, double t13sqsigmap = NAN, double t13sqsigmam = NAN)
		{
			t13sqbest = t13sqb;
			if (!std::isnan(t13sqsigmap)) t13sqsigmaplus = t13sqsigmap;
			if (!std::isnan(t13sqsigmam)) t13sqsigmaminus = t13sqsigmam;
		}
		void set23(double t23sqb, double t23sqsigmap = NAN, double t23sqsigmam = NAN)
		{
			t23sqbest = t23sqb;
			if (!std::isnan(t23sqsigmap)) t23sqsigmaplus = t23sqsigmap;
			if (!std::isnan(t23sqsigmam)) t23sqsigmaminus = t23sqsigmam;
		}
		void setdm21(double dm21sqb, double dm21sqsigmap = NAN, double dm21sqsigmam = NAN)
		{
			dm21sqbest = dm21sqb;
			if (!std::isnan(dm21sqsigmap)) dm21sqsigmaplus = dm21sqsigmap;
			if (!std::isnan(dm21sqsigmam)) dm21sqsigmaminus = dm21sqsigmam;
		}
		void setdm31(double dm31sqb, double dm31sqsigmap = NAN, double dm31sqsigmam = NAN)
		{
			dm31sqbest = dm31sqb;
			if (!std::isnan(dm31sqsigmap)) dm31sqsigmaplus = dm31sqsigmap;
			if (!std::isnan(dm31sqsigmam)) dm31sqsigmaminus = dm31sqsigmam;
		}
		void setdcp(double dcpb, double dcpsigmap = NAN, double dcpsigmam = NAN)
		{
			dcpbest = dcpb;
			if (!std::isnan(dcpsigmap)) dcpsigmaplus = dcpsigmap;
			if (!std::isnan(dcpsigmam)) dcpsigmaminus = dcpsigmam;
		}

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


	private:
		double t12sqbest, t13sqbest, t23sqbest, dm21sqbest, dm31sqbest, dcpbest;
		double t12sqsigmaplus, t13sqsigmaplus, t23sqsigmaplus, dm21sqsigmaplus, dm31sqsigmaplus, dcpsigmaplus;
		double t12sqsigmaminus, t13sqsigmaminus, t23sqsigmaminus, dm21sqsigmaminus, dm31sqsigmaminus, dcpsigmaminus;
		std::string ordering;
};

class oscillationexperiment
{
	public:
		oscillationexperiment(oscillationparams &osc){}
		~oscillationexperiment(){};
		virtual void setosc(oscillationparams &osc){}; //set oscillation parameters in an experiment
};

/* /////////////////////////////////////////////////////////////////////////////

DATA FROM INDIVIDUAL EXPERIMENTS BELOW

 ///////////////////////////////////////////////////////////////////////////// */

/* Data from the JUNO experiment. Numbers taken from xxxx.xxxxx p.x Tab. y  */
class JUNO : oscillationexperiment
{
	public:
		JUNO(oscillationparams &osc) : oscillationexperiment(osc)
		{
			setosc(osc);
		}
		~JUNO(){};
		void setosc(oscillationparams &osc)
		{
			//set the standard deviation to be the smaller one between the curret value and
			//the expected experimental value to be considered
			if (t12sqsigma < osc.get12plus()) osc.set12(osc.get12best(), t12sqsigma, osc.get12minus());
			if (t12sqsigma < osc.get12minus()) osc.set12(osc.get12best(), osc.get12plus(), t12sqsigma);
		}
	private:
		double t12sqsigma = 0.002;
};

int main()
{
	oscillationparams osc;
	std::cout << "before JUNO: sigma = " << osc.get12plus() << std::endl;
	JUNO JUNOEXP(osc);
	std::cout << "after JUNO: sigma = " << osc.get12plus() << std::endl;
	std::cout << "mass ordering: " << osc.getordering() << std::endl;
}
