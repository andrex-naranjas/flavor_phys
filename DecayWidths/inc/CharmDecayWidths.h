//CharmDecayWidths includes
#ifndef CHARMDECAYWIDTHS_H
#define CHARMDECAYWIDTHS_H

#include "ConfigSettings.h"

#include <string>
#include <vector>


class CharmDecayWidths{

public:
  CharmDecayWidths();
  virtual ~CharmDecayWidths();
  
  virtual void  initialize(Config config);
  virtual void  execute(Config config);

private:

  double MA; double MB; double MC;
  std::string modeExcitation="";
  double pi_val = 3.1415926536;
  double JA = 0.5;   std::vector<double> mJA;
  double LA = 1.;    std::vector<double> mLA;
  double SA = 0.5;   std::vector<double> mSA;
  double SB = 0.5;   std::vector<double> mSB;
  double SC = 0;     std::vector<double> mSC;
  double slight = 0; std::vector<double> m23; 
  double slightf = 1;std::vector<double> m24; 
  double s = 1.0;    std::vector<double> m;   
  double s1 = 0.5;   std::vector<double> m1;  
  double s2 = 0.5;   std::vector<double> m2;  
  double s3 = 0.5;   std::vector<double> m3;  
  double s4 = 0.5;   std::vector<double> m4;  
  double s5 = 0.5;   std::vector<double> m5;  

  virtual int KroneckerDelta(float i, float j);
  virtual std::vector<double> getMomentumProjections(double j_angular);
  virtual double ANGULAR_SUM(double alpha_d, double alpha_rho, double alpha_lam,
			     double alpha_mes, double k_value);
  virtual double DecayWidth(double decay, double gamma, double fi2_value, double angular_sum_value);

  virtual void PrintResults(Config config, std::vector<double> widths);

  virtual double EWCC(double MA, double MB, double MC);
  virtual double EB(double MA, double MB, double MC);
  virtual double K(double EB, double MB);
  virtual double FI2(double EB, double EWCC, double MA, double k_value);
  virtual double CBARFIN(double alpha_rho, double alpha_lam);
  virtual double CBARIN(double alpha_rho, double alpha_lam);
  virtual double CMESON(double alpha_mes);  
  virtual double C0(double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double F00(double alpha_d, double alpha_rho, double alpha_lam,double alpha_mes, double k_value);
  virtual double ARO0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double ALAM0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double BLAM0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double CLAM0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double A0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes);
  virtual double B0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double F0TOT(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value);

  virtual double I010(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I020(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value);


  //rho mode
  virtual double ARO0_rho(double alpha_lam, double alpha_mes);
  virtual double F01(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double BRO1(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double ARO1(double alpha_lam, double alpha_mes);
  virtual double ARO2(double alhpha_rho, double alpha_lam, double alpha_mes);
  virtual double BRO2(double alhpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double F0TOT_rho(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I01B0TOT(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  virtual double I02B0TOT(double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  

};

#endif //> !CHARMDECAYWIDTHS_H

// virtual double I01(double alpha_d, double alpha_rho, double alpha_lam,
// 		     double alpha_mes, double k_value, float i, float j);
// virtual double I20I10(double alpha_d, double alpha_rho, double alpha_lam,
// 			double alpha_mes, double k_value, float i, float j);

  // virtual double I01B0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  // virtual double I02B0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  // virtual double I03B0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  // virtual double I01B0TOT(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  // virtual double I02B0TOT(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
  // virtual double I03B0TOT(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value);
