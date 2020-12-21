//CharmDecay Class includes
#ifndef CHARMDECAYWIDTHS_CXX
#define CHARMDECAYWIDTHS_CXX

#include "CharmDecayWidths.h"
#include "WignerSymbols.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include <stdio.h>
#include <complex>


CharmDecayWidths::CharmDecayWidths()
{
}

CharmDecayWidths::~CharmDecayWidths(){}

void CharmDecayWidths::initialize(Config config)
{
  return;
}

void CharmDecayWidths::execute(Config config){

  std::vector<double> widths_results; widths_results.clear();//to store results
    
  int nStates = config.JA_quantum.size();

  for(int iStates = 0; iStates<nStates; iStates++){

    modeExcitation = config.ModeExc.at(iStates);
    //fetch quantum numbers and projections
    JA = config.JA_quantum.at(iStates);   mJA = getMomentumProjections(JA);
    LA = config.LA_quantum.at(iStates);   mLA = getMomentumProjections(LA);
    SA = config.SA_quantum.at(iStates);   mSA = getMomentumProjections(SA);
    SB = config.SB_quantum.at(iStates);   mSB = getMomentumProjections(SB);
    SC = config.SC_quantum.at(iStates);   mSC = getMomentumProjections(SC);
    s  = config.S_quantum.at(iStates);    m   = getMomentumProjections(s);
    s1 = config.S1_quantum.at(iStates);   m1  = getMomentumProjections(s1);
    s2 = config.S2_quantum.at(iStates);   m2  = getMomentumProjections(s2);
    s3 = config.S3_quantum.at(iStates);   m3  = getMomentumProjections(s3);
    s4 = config.S4_quantum.at(iStates);   m4  = getMomentumProjections(s4);
    s5 = config.S5_quantum.at(iStates);   m5  = getMomentumProjections(s5);
    slight  = config.Slight_quantum.at(iStates); m23 = getMomentumProjections(slight);
    slightf = config.Slightf_quantum.at(iStates);m24 = getMomentumProjections(slightf);
    
    MA = config.MassA.at(iStates);
    MB = config.MassB.at(iStates);
    MC = config.MassC.at(iStates);

    double EB_value = EB(MA,MB,MC);
    
    double gamma     = 9.22051;
    double alpha_d   = 0.0;
    double alpha_rho = 0.437553;
    double alpha_lam = 0.523495;
    double alpha_mes = 0.0;
    double decay_coef= 0.0;
    if(config.decayProd.at(iStates) == "1st")      {alpha_mes = 0.40; decay_coef = 0.125;}
    else if(config.decayProd.at(iStates) == "2nd") {alpha_mes = 0.46; decay_coef = 0.250;}
      
    double k_value; k_value = K(EB_value, MB);
    double EWCC_value = EWCC(MA, MB, MC);
    
    double sum_value = ANGULAR_SUM(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value);
    double fi2_value = FI2(EB_value, EWCC_value, MA, k_value);
    double decayWidth = DecayWidth(decay_coef, gamma, fi2_value, sum_value);
    widths_results.push_back(decayWidth);
  }
  if(config.PrintResults=="True")
    PrintResults(config, widths_results);
  return;
}

double CharmDecayWidths::DecayWidth(double dec_coef, double gamma, double fi2_value, double angular_sum_value){
  double GeV = 1000.;
  double decayWidth = dec_coef * std::pow(gamma, 2) * fi2_value * (1./(2*JA + 1)) * angular_sum_value;
  return decayWidth*GeV;
}

double CharmDecayWidths::ANGULAR_SUM(double alpha_d, double alpha_rho, double alpha_lam,
				     double alpha_mes, double k_value){
  
  WignerSymbols *m_wigner = new WignerSymbols();
  double outerSum = 0;
  double finalIntegral1=0., finalIntegral2=0.;

  if(modeExcitation == "lam"){
    finalIntegral1=I010(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I020(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value);
  }else if(modeExcitation == "rho"){
    finalIntegral1=I01B0TOT(alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I02B0TOT(alpha_rho, alpha_lam, alpha_mes, k_value);
  }
  
  for(int iMJA = 0; iMJA<mJA.size(); iMJA++){
    
    double innerSum = 0;    
    for(int iMLA = 0; iMLA<mLA.size(); iMLA++)
      for(int iMSA = 0; iMSA<mSA.size(); iMSA++)
	for(int iM = 0; iM<m.size(); iM++)
	  for(int iM24 = 0; iM24<m24.size(); iM24++)
	    for(int iM1 = 0; iM1<m1.size(); iM1++)
	      for(int iMSB = 0; iMSB<mSB.size(); iMSB++)
		for(int iM3 = 0; iM3<m3.size(); iM3++)
		  for(int iM5 = 0; iM5<m5.size(); iM5++)
		    for(int iMSC = 0; iMSC<mSC.size(); iMSC++)
		      for(int iM23 = 0; iM23<m23.size(); iM23++)
			for(int iM4 = 0; iM4<m4.size(); iM4++)
			  for(int iM2 = 0; iM2<m2.size(); iM2++){
			    int delta1=0, delta2=0;
			    if(modeExcitation == "lam"){
			      delta1 = KroneckerDelta(m.at(iM), mLA.at(iMLA));
			      delta2 = KroneckerDelta(m.at(iM), 0)*KroneckerDelta(mLA.at(iMLA), 0);
			    }else if(modeExcitation == "rho"){
			      delta1 = KroneckerDelta(m.at(iM), 0)*KroneckerDelta(mLA.at(iMLA), 0);
			      delta2 = KroneckerDelta(m.at(iM), mLA.at(iMLA));}
			    
			    double dummy = (finalIntegral1*delta1 + finalIntegral2*delta2)*
			      m_wigner->wigner3j(LA, SA, JA, mLA.at(iMLA), mSA.at(iMSA), (-1.0)*mJA.at(iMJA))*
			      m_wigner->wigner3j(1, 1, 0, m.at(iM), (-1.0)*m.at(iM), 0)*
			      m_wigner->wigner3j(slightf, 0.5, SB, m24.at(iM24), m1.at(iM1), (-1.0)*mSB.at(iMSB))*
			      m_wigner->wigner3j(0.5, 0.5, SC, m3.at(iM3), m5.at(iM5), (-1.0)*mSC.at(iMSC))*
			      m_wigner->wigner3j(slight, 0.5, SA, m23.at(iM23), m1.at(iM1), (-1.0)*mSA.at(iMSA))*
			      m_wigner->wigner3j(0.5, 0.5, 1, m4.at(iM4), m5.at(iM5), m.at(iM))*
			      m_wigner->wigner3j(0.5, 0.5, slightf, m2.at(iM2), m4.at(iM4), (-1.0)*m24.at(iM24))*
			      m_wigner->wigner3j(0.5,0.5,slight, m2.at(iM2), m3.at(iM3), (-1.0)*m23.at(iM23))*
			      std::pow(3 * (2*JA+1) * (2*slight+1) * (2*slightf+1) * (2*SA+1) * (2*SB+1) * (2*SC+1),0.5)*
			      std::pow(-1.0,SA-LA-mJA.at(iMJA))*
			      std::pow(-1.0,1+m.at(iM)-mSA.at(iMSA)-mSB.at(iMSB)-slight-slightf-m23.at(iM23)-m24.at(iM24)-mSC.at(iMSC));			      
			      innerSum+=dummy;			      
			      dummy = 0;
			      }
    outerSum += std::pow(innerSum,2);        
  }
  return outerSum;
}

std::vector<double> CharmDecayWidths::getMomentumProjections(double j_angular){
  //gets the m projections "m_projection" for a given angular momentum "j_angular"
  std::vector<double> angularProjections; angularProjections.clear();
  if(j_angular==0){ angularProjections.push_back(0); return angularProjections;}
  
  double m_projection = (-1.0)*j_angular;
  do{
    angularProjections.push_back(m_projection);
    m_projection++;        
  }while(m_projection<=j_angular);
  
  return angularProjections;
}

int CharmDecayWidths::KroneckerDelta(float i, float j){
  if(i==j) return 1;
  else return 0;
}

double CharmDecayWidths::EWCC(double MA, double MB, double MC){
  double value = (0.5*(std::pow(MA,2) + std::pow(MC,2) - std::pow(MB,2)) )/ MA;
  return value;
}

double CharmDecayWidths::EB(double MA, double MB, double MC){  
  double value = ( 0.5*(std::pow(MA,2) - std::pow(MC,2) + std::pow(MB,2)) )/ MA;
  return value;
}

double CharmDecayWidths::K(double EB, double MB){
  double value = std::sqrt(std::pow(EB,2) - std::pow(MB,2) );
  return value;
}

double CharmDecayWidths::FI2(double EB, double EWCC, double MA, double k_value){
  double value = (2*(pi_val)*k_value*(EB * EWCC)) / MA ; 
  return value;
}

double CharmDecayWidths::CBARFIN(double alpha_rho, double alpha_lam){
  double value1 = std::pow(3,0.75) * std::pow( 8.0/(3*std::pow(pi_val, 0.5)), 0.5);
  double value2 = std::pow(1.0/( std::pow(alpha_lam,2) ), 1.25);
  double value3 = std::pow(1.0/( pi_val * std::pow(alpha_rho, 2) ), 0.75);
  return value1 * value2 * value3;// * mycomplex;//define complex if needed
}

double CharmDecayWidths::CBARIN(double alpha_rho, double alpha_lam){
  double value1 = std::pow(3,0.75);
  double value2 = std::pow(1.0/( pi_val * std::pow(alpha_rho, 2) ), 0.75);
  double value3 = std::pow(1.0/( pi_val * std::pow(alpha_lam, 2) ), 0.75);
  return value1 * value2 * value3 ;
}

double CharmDecayWidths::CMESON(double alpha_mes){
  double value = std::pow(1.0/( pi_val * std::pow(alpha_mes, 2) ), 0.75);
  return value;
}

double CharmDecayWidths::C0(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1.0 / (3*std::pow(3,0.5));
  double value2 = CBARFIN(alpha_rho, alpha_lam) * CBARIN(alpha_rho, alpha_lam) * CMESON(alpha_mes);
  return value1 * value2;
}

double CharmDecayWidths::F00(double alpha_d, double alpha_rho, double alpha_lam,
			     double alpha_mes, double k_value){
  double value1 = k_value * (-6.0 * std::pow(alpha_d, 2) );
  double value2 = std::pow(k_value,2)*( 4.0/std::pow(alpha_lam,2) + 3.0/std::pow(alpha_mes,2) + 4.0/std::pow(alpha_rho,2));
  return (1./24.)*(value1+value2);  
}

double CharmDecayWidths::ARO0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = (1.0/6.0) * (1./std::pow(alpha_lam,2));
  double value2 = (1.0/4.0) * (1./std::pow(alpha_mes,2));
  double value3 = (5.0/6.0) * (1./std::pow(alpha_rho,2));
  return std::pow(value1+value2+value3,0.5);
}

double CharmDecayWidths::ALAM0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = (1./12.)*( 12./ std::pow(alpha_lam,2) + 1./std::pow(alpha_mes,2) );
  double value2 = 1./( 48 * std::pow(alpha_mes,4) * std::pow(ARO0(alpha_d, alpha_rho, alpha_lam, alpha_mes), 2) );
  return std::pow(value1-value2, 0.5);
}

double CharmDecayWidths::BLAM0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1 = std::pow(alpha_d,2) / (4 * std::pow(6,0.5));
  double value2 = k_value / ( std::pow(6,0.5) * std::pow(alpha_lam,2) );
  double value3 = k_value / (2 * std::pow(6,0.5) * std::pow(alpha_mes,2) );
  double value4p= k_value * ( -2.0/std::pow(alpha_lam,2) - 3.0/std::pow(alpha_mes,2) - 4.0/std::pow(alpha_rho,2));
  double value4pp= 24 * std::pow(6,0.5) * std::pow(alpha_mes,2) * std::pow(ARO0(alpha_d, alpha_rho, alpha_lam, alpha_mes), 2) ;
  double value4 = value4p/value4pp;
  return (value1 - value2 - value3 - value4)/(2.0 * ALAM0(alpha_d, alpha_rho, alpha_lam, alpha_mes) );
}

double CharmDecayWidths::CLAM0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value){

  double value1p= std::pow(ARO0(alpha_d, alpha_rho, alpha_lam, alpha_mes), 2) ;
  double value1 = std::pow(alpha_d,4)/ (128 * value1p);
  double value2p= -2.0/std::pow(alpha_lam,2) - 3.0/std::pow(alpha_mes,2) - 4.0/std::pow(alpha_rho,2);
  double value2 = (k_value * std::pow(alpha_d,2) * value2p) / (96 * value1p);
  double value3 = (std::pow(k_value,2) * std::pow(value2p,2))/(288 * value1p);
  return -value1-value2-value3;
}

double CharmDecayWidths::A0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1./(2 * std::pow(3,0.5) * std::pow(alpha_mes,2));
  double value2 = 2. * ARO0(alpha_d, alpha_rho, alpha_lam, alpha_mes) ;
  return value1/value2;
}

double CharmDecayWidths::B0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::pow(alpha_d,2)/(4* std::pow(2,0.5));
  double value2 = k_value /( 3 * std::pow(2,0.5) * std::pow(alpha_lam,2) );
  double value3 = k_value /( 2 * std::pow(2,0.5) * std::pow(alpha_mes,2) );
  double value4 = (k_value * std::pow(2,0.5)) / ( 3 * std::pow(alpha_rho,2) );
  double value5 = 2. * ARO0(alpha_d, alpha_rho, alpha_lam, alpha_mes) ;
  return (value1-value2-value3-value4)/value5;
}

double CharmDecayWidths::F0TOT(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F00(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2 = std::pow(BLAM0(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value),2);
  double value3 = CLAM0(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1-value2+value3;
}

//we start here
double CharmDecayWidths::I010(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0(alpha_rho, alpha_lam, alpha_mes);
  double value2 = std::exp((-1.0)*F0TOT(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value));
  double value3 = std::pow(( std::pow(pi_val,0.5) / ARO0(alpha_d, alpha_rho, alpha_lam, alpha_mes )), 3);
  double value4 = (3 * std::pow(pi_val,0.5))/(8 * std::pow(ALAM0(alpha_d, alpha_rho, alpha_lam, alpha_mes),5));
  double value5 = std::pow(6,0.5)/3. - (std::pow(2,0.5) * A0(alpha_d, alpha_rho, alpha_lam, alpha_mes))/ (ARO0(alpha_d, alpha_rho, alpha_lam, alpha_mes));
  return value1*value2*value3*(value4*(value5));
}

double CharmDecayWidths::I020(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1 = C0(alpha_rho, alpha_lam, alpha_mes);
  double value2 = A0(alpha_d,  alpha_rho,  alpha_lam,  alpha_mes);
  double value3 = B0(alpha_d,  alpha_rho,  alpha_lam,  alpha_mes,  k_value);
  double value4 = BLAM0(alpha_d, alpha_rho,alpha_lam,  alpha_mes,  k_value);
  double value5 = ALAM0(alpha_d, alpha_rho,alpha_lam,  alpha_mes);
  double value6 = ARO0(alpha_d,  alpha_rho,  alpha_lam,  alpha_mes);
  double value7 = std::exp((-1.0)*F0TOT(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value));

  double valueA = std::pow( std::pow(pi_val,0.5) / value6, 3);
  double valueB = (0.75*std::pow(pi_val,0.5))*(value4/value5);
  double valueC = 1./std::pow(value5,3);
  double valueD = std::pow(2,0.5)*( value3/value6 - (value2*value4)/(value6*value5) );
  double valueE = (1./3.)*std::pow(6,0.5)*(value4/value5) + 2*k_value;

  return value1 * value7 * valueA * valueB * valueC * (valueD + valueE); 
}


// for rho harmonic oscilator modes
double CharmDecayWidths::ARO0_rho(double alpha_lam, double alpha_mes){
  double value1 = 1./std::pow(alpha_lam,2);
  double value2 = 1./(12*std::pow(alpha_mes, 2));  
  return std::pow(value1+value2,0.5);
}

double CharmDecayWidths::F01(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = 1./(2*std::pow(6,0.5)*alpha_lam*alpha_rho) + 1./(4*std::pow(6,0.5)*alpha_mes*alpha_mes);  
  double value1pp = value1p*value1p;
  double value1ppp= 1./std::pow(ARO0_rho(alpha_lam, alpha_mes),2);
  double value1 = value1pp*value1ppp;
  double value2 = 1./(8*std::pow(alpha_mes,2));  
  return k_value*k_value*(value1-value2);
}

double CharmDecayWidths::BRO1(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = 1./(2*alpha_lam*alpha_rho) + 1./(4*alpha_mes*alpha_mes);
  double value1pp = 1./(3*std::pow(ARO0_rho(alpha_lam, alpha_mes), 2));
  double value1 = value1p * value1pp;
  double value2 = 2.;
  return k_value*(value1-value2);
}

double CharmDecayWidths::ARO1(double alpha_lam, double alpha_mes){
  double value1 = std::pow(2,0.5);
  double value2 = 1./(6 * std::pow(2,0.5) * std::pow(ARO0_rho(alpha_lam, alpha_mes),2) * std::pow(alpha_mes,2) );
  return -value1+value2;
}

double CharmDecayWidths::ARO2(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1p  = 4./std::pow(alpha_rho,2);
  double value1pp = 1./std::pow(alpha_mes,2);
  double value1   = 0.25*(value1p + value1pp);
  double value2   = 1./(48*std::pow(alpha_mes,4)*std::pow(ARO0_rho(alpha_lam, alpha_mes),2));
  return std::pow(value1-value2, 0.5);
}

double CharmDecayWidths::BRO2(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1    = k_value/(4*std::pow(2,0.5)*std::pow(alpha_mes,2));
  double value2    = k_value/(2*std::pow(2,0.5)*alpha_rho*alpha_lam);
  double value3p   = k_value/(4*std::pow(3,0.5)*std::pow(alpha_mes,2));
  double value3ppp =      1./(2*std::pow(6,0.5)*alpha_rho*alpha_lam) + 1./(4*std::pow(6,0.5)*std::pow(alpha_mes,2));
  double value3pp  = 1./std::pow(ARO0_rho(alpha_lam, alpha_mes),2);
  double value3 = value3p * value3pp * value3ppp;
  return -value1 - value2 + value3;
}

double CharmDecayWidths::F0TOT_rho(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2 = std::pow(BRO2(alpha_rho, alpha_lam, alpha_mes, k_value)/ARO2(alpha_rho, alpha_lam, alpha_mes), 2);
  double value3 = (1./3.)*(std::pow(k_value,2)*(1./(alpha_rho*alpha_lam)));
  return value1 + value2 - value3;
}

double CharmDecayWidths::I01B0(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = C0(alpha_rho, alpha_lam, alpha_mes) / std::pow(ARO0_rho(alpha_lam, alpha_mes), 3) ;
  double value1pp = BRO2(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2(alpha_rho, alpha_lam, alpha_mes), 2);
  double value1 = value1p*value1pp;
  double value2 = value1pp*ARO1(alpha_lam, alpha_mes) + BRO1(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value3 = 1./(4.* std::pow( ARO2(alpha_rho, alpha_lam, alpha_mes), 3));
  return value1*value2*value3*pi_val*pi_val;
}

double CharmDecayWidths::I02B0(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0(alpha_rho, alpha_lam, alpha_mes) / std::pow(ARO0_rho(alpha_lam, alpha_mes), 3);
  double value2 = (3. * ARO1(alpha_lam, alpha_mes) ) / (8.* std::pow(  ARO2(alpha_rho, alpha_lam, alpha_mes), 5) );
  return pi_val*pi_val*value1*value2;
}

double CharmDecayWidths::I01B0TOT(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_rho(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double CharmDecayWidths::I02B0TOT(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_rho(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

void CharmDecayWidths::PrintResults(Config config, std::vector<double> widths){  
  std::string exc=""; double decay_val=0;
  printf("-----------------------------------------------------------------------------------------\n");
  printf("   MA       MB        MC       ");
  printf("JA     LA     SA     SB     Slight    ");
  printf("Exct.    DecayWidth \n");
  printf("-----------------------------------------------------------------------------------------\n");
  
  int nStates = config.JA_quantum.size();
  for(int iStates = 0; iStates<nStates; iStates++){    
    JA = config.JA_quantum.at(iStates);
    LA = config.LA_quantum.at(iStates);
    SA = config.SA_quantum.at(iStates);
    SB = config.SB_quantum.at(iStates);
    SC = config.SC_quantum.at(iStates);
    slight  = config.Slight_quantum.at(iStates); 
    slightf = config.Slightf_quantum.at(iStates);    
    MA = config.MassA.at(iStates);
    MB = config.MassB.at(iStates);
    MC = config.MassC.at(iStates);
    exc = config.ModeExc.at(iStates);
    decay_val = widths.at(iStates);    
    printf("%.4f    %.4f    %.4f    %.1f    %.1f    %.1f     %.1f      %.1f     %s      %.6f\n",
	     MA,      MB,     MC,      JA,    LA,    SA,      SB,      slight,  exc.c_str(),  decay_val);
  }
  printf("-----------------------------------------------------------------------------------------\n");
  return;
}


#endif

// double CharmDecayWidths::I20I10(double alpha_d, double alpha_rho, double alpha_lam,
// 				double alpha_mes, double k_value, float i, float j){
//   return I020(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value)*KroneckerDelta(i,0)*KroneckerDelta(j,0);
// }

// double CharmDecayWidths::I01(double alpha_d, double alpha_rho, double alpha_lam,
// 			     double alpha_mes, double k_value, float i, float j){
//   return I010(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value)*KroneckerDelta(i,j);
// }

// double CharmDecayWidths::BRO0(double alpha_d, double alpha_rho, double alpha_lam,
// 			      double alpha_mes, double k_value){
//   double p_lam=0;
  
//   double value1 = std::pow(alpha_d,2) / (4 * std::pow(2,0.5));
//   double value2 = k_value / (3 * std::pow(2,0.5) * std::pow(alpha_lam,2) );
//   double value3 = k_value / (2 * std::pow(2,0.5) * std::pow(alpha_mes,2) );
//   double value4 = p_lam   / (2 * std::pow(3,0.5) * std::pow(alpha_mes,2) );
//   double value5 = std::pow(2,0.5) * k_value / (3 * std::pow(alpha_rho,2) );
//   return (0.5*(value1-value2-value3+value4-value4))/ARO0(alpha_d, alpha_rho, alpha_lam, alpha_mes);
// }

// double CharmDecayWidths::CRO0(double alpha_d, double alpha_rho, double alpha_lam,
// 			      double alpha_mes, double k_value){
//   double p_lam=0;

//   double value1 = (p_lam * std::pow(alpha_d,2)) / (4 * std::pow(6,0.5));
//   double value2 = (k_value * p_lam) / (std::pow(6,0.5) * std::pow(alpha_lam, 2) );
//   double value3 = (1./24.)*std::pow(p_lam,2)*( 24./ std::pow(alpha_lam,2) + 2./std::pow(alpha_mes,2) );
//   double value4 = (k_value * p_lam) / (2 * std::pow(6,0.5) * std::pow(alpha_mes, 2) );
//   return value1 - value2 + value3  - value4;
// }

// double CharmDecayWidths::I01B0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
//   double value1 = C0(alpha_rho, alpha_lam, alpha_mes);
//   double value2 = (3 * std::pow(pi_val,0.5))/(8 * std::pow(ALAM0(alpha_d, alpha_rho, alpha_lam, alpha_mes),5));
//   double value3 = std::pow(std::pow(pi_val,0.5) / ARO0(alpha_d, alpha_rho, alpha_lam, alpha_mes),3);
//   double value4 = (std::pow(6,0.5) * A0(alpha_d, alpha_rho, alpha_lam, alpha_mes))/ (3 * ARO0(alpha_d, alpha_rho, alpha_lam, alpha_mes));
//   double value5 = std::pow(2,0.5)*std::pow(A0(alpha_d, alpha_rho, alpha_lam, alpha_mes)/ARO0(alpha_d, alpha_rho, alpha_lam, alpha_mes),2);
//   return value1*value2*(value5-value4);
// }

// double CharmDecayWidths::I02B0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
//   double value1 = C0(alpha_rho, alpha_lam, alpha_mes);
//   double value2 = A0(alpha_d,  alpha_rho,  alpha_lam,  alpha_mes);
//   double value3 = B0(alpha_d,  alpha_rho,  alpha_lam,  alpha_mes,  k_value);
//   double value4 = BLAM0(alpha_d, alpha_rho,alpha_lam,  alpha_mes,  k_value);
//   double value5 = ALAM0(alpha_d, alpha_rho,alpha_lam,  alpha_mes);
//   double value6 = ARO0(alpha_d,  alpha_rho,  alpha_lam,  alpha_mes);
//   double valueA = (std::pow(6,0.5)/3.0 ) * (value4/value5);
//   double valueB = std::pow(2,0.5) * ( (value2*value4) / (value6*value5 ) );
//   double valueC = std::pow(2,0.5) * (value3 / value6 ) + 2*k_value;
//   double valueD = (-1.0)*(value2*value4 / value6*value5 );
//   double valueE = value3/value6;
//   double valueG =  (3*pi_val*pi_val/4)*std::pow(1./(value5*value6),2);

//   return value1 * (valueA-valueB+valueC)*(valueD+valueE)*valueG;
// }

// double CharmDecayWidths::I03B0(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
//    double value1 = C0(alpha_rho, alpha_lam, alpha_mes);
//    double value2 = std::pow(ARO0(alpha_d,  alpha_rho,  alpha_lam,  alpha_mes),5);
//    double value3 = std::pow(ALAM0(alpha_d, alpha_rho,alpha_lam,  alpha_mes),3);
//    return (value1 * pi_val * pi_val * (3./8.))/(value2 * value3) ;
// }


// double CharmDecayWidths::I01B0TOT(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
//    double value1 = I01B0(alpha_d, alpha_rho, alpha_lam,  alpha_mes,  k_value);
//    double value2 = std::exp(F0TOT( alpha_d,  alpha_rho,  alpha_lam,  alpha_mes,  k_value));
//    return value1*value2;
// }

// double CharmDecayWidths::I02B0TOT(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
//    double value1 = I02B0(alpha_d, alpha_rho, alpha_lam,  alpha_mes,  k_value);
//    double value2 = std::exp(F0TOT( alpha_d,  alpha_rho,  alpha_lam,  alpha_mes,  k_value));
//    return value1*value2;
// }

// double CharmDecayWidths::I03B0TOT(double alpha_d, double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
//    double value1 = I03B0(alpha_d, alpha_rho, alpha_lam,  alpha_mes,  k_value);
//    double value2 = std::exp(F0TOT( alpha_d,  alpha_rho,  alpha_lam,  alpha_mes,  k_value));
//    return value1*value2;
// }
