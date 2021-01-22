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

double CharmDecayWidths::execute(double ma_val, double mb_val, double mc_val, double sa_val,
				 double la_val, double ja_val, double sl_val,
				 int baryon, int excMode, int prodDecay){  
  // decay product masses
  MA = ma_val;
  MB = mb_val;
  MC = mc_val;
  if(MA<MB+MC) return 0; //energy conservation

  // which baryon, mode, decay product
  int baryonFlag = baryon;
  modeExcitation = excMode;
  int decayProd  = prodDecay;
  double alpha_rho = 0.,alpha_lam = 0.,alpha_mes = 0.,flav_coup= 0.;
  double slf_val=0., sb_val=0.;

  // options according the type of baryon
  if(baryonFlag==1)      {alpha_rho = 0.458724; alpha_lam = 0.540131;} //omegas
  else if(baryonFlag==2) {alpha_rho = 0.437553; alpha_lam = 0.523495;} //cacades 6-plet
  else if(baryonFlag==3) {alpha_rho = 0.412767; alpha_lam = 0.502336;} //sigmas
  else if(baryonFlag==4) {alpha_rho = 0.412767; alpha_lam = 0.502336;} //lambdas
  else if(baryonFlag==5) {alpha_rho = 0.437553; alpha_lam = 0.523495;} //cascades 3-plet
  
  // options according the decay products
  if(baryonFlag==1){// omegas
    if(decayProd == 1)    {slf_val=0.; sb_val=0.5; flav_coup = 1./3.;} //Xi+K
    else if(decayProd==2) {slf_val=1.; sb_val=0.5; flav_coup = 1./3.;} //Xi'+K
    else if(decayProd==3) {slf_val=1.; sb_val=1.5; flav_coup = 1./3.;} //Xi*+K
  }else if(baryonFlag==2 or baryonFlag==5){//cascades (6-plet, 3-plet)
    if(decayProd == 1)    {slf_val=0.; sb_val=0.5; flav_coup = 1./12.;}//lambda+K
    else if(decayProd==2) {slf_val=0.; sb_val=0.5; flav_coup = 1./8.;} //Xi+Pi
    else if(decayProd==3) {slf_val=1.; sb_val=0.5; flav_coup = 1./8.;} //Xi'(2578)+Pi comp. w/mathematica
    else if(decayProd==4) {slf_val=1.; sb_val=1.5; flav_coup = 1./8.;} //Xi*+pi
    else if(decayProd==5) {slf_val=1.; sb_val=0.5; flav_coup = 1./4.;} //Sigma+K comp. w/ mathematica
    else if(decayProd==6) {slf_val=1.; sb_val=1.5; flav_coup = 1./4.;} //Sigma*+K
    else if(decayProd==7) {slf_val=0.; sb_val=0.5; flav_coup = 1./36.;}//Xi+eta
  }else if(baryonFlag==3){// sigmas
    if(decayProd == 1)    {slf_val=1.; sb_val=0.5; flav_coup = 1./3.;} //Sigma+Pi
    else if(decayProd==2) {slf_val=1.; sb_val=1.5; flav_coup = 1./3.;} //Sigma*+pi
    else if(decayProd==3) {slf_val=0.; sb_val=0.5; flav_coup = 1./2.;} //Lambda+Pi
    else if(decayProd==4) {slf_val=1.; sb_val=0.5; flav_coup = 1./18.;}//Sigma+eta
    else if(decayProd==5) {slf_val=0.; sb_val=0.5; flav_coup = 3./2.;} //Xi+K
  }else if(baryonFlag==4){// Lambdas
    if(decayProd == 1)    {slf_val=1.; sb_val=0.5; flav_coup = 1./2.;} //Sigma+Pi
    else if(decayProd==2) {slf_val=1.; sb_val=1.5; flav_coup = 1./2.;} //Sigma*+pi
    else if(decayProd==3) {slf_val=0.; sb_val=0.5; flav_coup = 1./18.;}//Lambda+eta
  }

  alpha_mes = ALPHA_MES(MC);

  //test  
  // slf_val=1.; sb_val=0.5; flav_coup = 1./4.;
  // MA = 3.080;
  // MB = 2.286;
  // MC = 0.493;
  // alpha_rho = 0.437553;
  // alpha_lam = 0.523495;
  // alpha_mes = 0.46;
  
  //fetch quantum numbers and projections
  SA = sa_val;      mSA = getMomentumProjections(SA); 
  LA = la_val;      mLA = getMomentumProjections(LA);
  JA = ja_val;      mJA = getMomentumProjections(JA);  
  SB = sb_val;      mSB = getMomentumProjections(SB);
  slight = sl_val;  m23 = getMomentumProjections(slight);
  slightf= slf_val; m24 = getMomentumProjections(slightf);

  //values are the same for all states (at least for now!)
  SC = 0.0;   mSC = getMomentumProjections(SC);
  s  = 1.0;   m   = getMomentumProjections(s);
  s1 = 0.5;   m1  = getMomentumProjections(s1);
  s2 = 0.5;   m2  = getMomentumProjections(s2);
  s3 = 0.5;   m3  = getMomentumProjections(s3);
  s4 = 0.5;   m4  = getMomentumProjections(s4);
  s5 = 0.5;   m5  = getMomentumProjections(s5);

  double EB_value = EB(MA,MB,MC);
  
  double gamma     = 17.25/(std::pow(7,0.5));
  double alpha_d   = 0.0;
  //gamma = 9.22051;
      
  double k_value; k_value = K(EB_value, MB);
  double EWCC_value = EWCC(MA, MB, MC);
  
  double sum_value  = ANGULAR_SUM(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value);
  double fi2_value  = FI2(EB_value, EWCC_value, MA, k_value);
  double decayWidth = DecayWidth(flav_coup, gamma, fi2_value, sum_value);  
      
  return decayWidth;
}

double CharmDecayWidths::DecayWidth(double flav_coup, double gamma, double fi2_value, double angular_sum_value){
  double GeV = 1000.;
  double decayWidth = flav_coup * std::pow(gamma, 2) * fi2_value * (1./(2*JA + 1)) * angular_sum_value;
  return decayWidth*GeV;
}

double CharmDecayWidths::ALPHA_MES(double MC_value){
  if(MC_value == 0.140) return 0.40;
  else return 0.46;
}

double CharmDecayWidths::ANGULAR_SUM(double alpha_d, double alpha_rho, double alpha_lam,
				     double alpha_mes, double k_value){
  
  WignerSymbols *m_wigner = new WignerSymbols();
  double outerSum = 0;
  double finalIntegral1=0., finalIntegral2=0.;

  if(modeExcitation == 1 && LA==1){//P-WAVE
    finalIntegral1=I010(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I020(alpha_d, alpha_rho, alpha_lam, alpha_mes, k_value);
  }else if(modeExcitation == 2 && LA==1){//P-WAVE
    finalIntegral1=I01B0TOT(alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I02B0TOT(alpha_rho, alpha_lam, alpha_mes, k_value);
  }else if(modeExcitation == 1 && LA==2){//D-WAVE
    finalIntegral1=I01B0TOT_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I02B0TOT_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  }else if(modeExcitation == 1 && LA==0){//RADIAL
    finalIntegral1=I01B0TOT_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
    finalIntegral2=I02B0TOT_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
  }
  
  for(int iMJA = 0; iMJA<(int)mJA.size(); iMJA++){
    
    double innerSum = 0;
    for(int iMLA = 0; iMLA<(int)mLA.size(); iMLA++)
      for(int iMSA = 0; iMSA<(int)mSA.size(); iMSA++)
	for(int iM = 0; iM<(int)m.size(); iM++)
	  for(int iM24 = 0; iM24<(int)m24.size(); iM24++)
	    for(int iM1 = 0; iM1<(int)m1.size(); iM1++)
	      for(int iMSB = 0; iMSB<(int)mSB.size(); iMSB++)
		for(int iM3 = 0; iM3<(int)m3.size(); iM3++)
		  for(int iM5 = 0; iM5<(int)m5.size(); iM5++)
		    for(int iMSC = 0; iMSC<(int)mSC.size(); iMSC++)
		      for(int iM23 = 0; iM23<(int)m23.size(); iM23++)
			for(int iM4 = 0; iM4<(int)m4.size(); iM4++)
			  for(int iM2 = 0; iM2<(int)m2.size(); iM2++){
			    int delta1=0, delta2=0;
			    if(modeExcitation == 1 && LA == 1){//P-WAVE
			      delta1 = KroneckerDelta(m.at(iM), mLA.at(iMLA));
			      delta2 = KroneckerDelta(m.at(iM), 0)*KroneckerDelta(mLA.at(iMLA), 0);
			    }else if(modeExcitation == 2 && LA ==1){//P-WAVE
			      delta1 = KroneckerDelta(m.at(iM), 0)*KroneckerDelta(mLA.at(iMLA), 0);
			      delta2 = KroneckerDelta(m.at(iM), mLA.at(iMLA));
			    }else if(modeExcitation == 1 && LA == 2){//D-WAVE		      
			      delta1 = KroneckerDelta(m.at(iM), 0) * KroneckerDelta(mLA.at(iMLA), 0);
			      delta2 = m_wigner->wigner3j(1, 1, LA, m.at(iM), 0, (-1.0)*m.at(iM)) *
				std::pow(-1.0, (-1.0)*m.at(iM)) * std::pow((2*LA+1),0.5) * KroneckerDelta(m.at(iM), mLA.at(iMLA));
			    }else if(modeExcitation == 1 && LA == 0){//RADIAL
			      delta1 = KroneckerDelta(m.at(iM), 0);
			      delta2 = KroneckerDelta(m.at(iM), 0);}
				
			    			    
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

double CharmDecayWidths::I02B0(double alpha_rho, double alpha_lam, double alpha_mes){
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
  double value2 = I02B0(alpha_rho, alpha_lam, alpha_mes);
  return value1*value2;
}


// D-wave excitations
// lambda wave functions
double CharmDecayWidths::CBARFIN_DWAVE(double alpha_rho, double alpha_lam){
  double value1 = (std::pow(3,0.75) * std::pow(16,0.5))/(std::pow(15,0.5) * (std::pow(1.772453667,0.5)) );
  double value2 = std::pow(( 1.0/alpha_lam ), 3.5);
  double value3 = std::pow(1.0/( pi_val * std::pow(alpha_rho, 2) ), 0.75);
  return value1 * value2 * value3;// * mycomplex;//define complex if needed
}

//CBARIN_DWAVE --> same as P-wave (CBARIN)
//CMESON_DWAVE --> same as P-wave (CMESON)

double CharmDecayWidths::C0_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1.0 / (3.0*std::pow(3,0.5));
  double value2 = CBARFIN_DWAVE(alpha_rho, alpha_lam) * CBARIN(alpha_rho, alpha_lam) * CMESON(alpha_mes);
  return value1 * value2;
}

double CharmDecayWidths::ARO0_DWAVE(double alpha_rho, double alpha_mes){
  double value1 = 1.0/std::pow(alpha_rho,2);  
  double value2 = 1.0/(4.0*std::pow(alpha_mes,2));
  return std::pow(value1+value2,0.5);
}

double CharmDecayWidths::F0_DWAVE(double k_value, double alpha_rho, double alpha_lam){
  double value1 = std::pow(k_value, 2) /(alpha_rho * alpha_lam) ;
  return (-1.0/3.0) * value1;
}

double CharmDecayWidths::F01_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  
  double value1p  = 1.0/(2.0*std::pow(2,0.5)*alpha_lam*alpha_rho) + 1./(4*std::pow(2,0.5)*alpha_mes*alpha_mes);
  
  double value1pp = value1p*value1p;
  
  double value1ppp= 1./(std::pow(ARO0_DWAVE(alpha_rho, alpha_mes),2));
  
  double value1 = value1pp*value1ppp;
  
  double value2 = 1.0/(8.0*std::pow(alpha_mes,2));
  
  return std::pow(k_value,2)*(value1-value2);
}

double CharmDecayWidths::BRO2_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = 1.0/(2.0*alpha_lam*alpha_rho) + 1.0/(4.0*alpha_mes*alpha_mes);
  double value2 = 1.0/(std::pow(ARO0_DWAVE(alpha_rho, alpha_mes),2));
  return k_value*(value1*value2 - 2.0);
}

double CharmDecayWidths::ARO1_DWAVE(double alpha_rho, double alpha_mes){
  double value1 = 0.816496581;
  double value2 = 1./(2.0 * std::pow(6,0.5) * std::pow(ARO0_DWAVE(alpha_rho, alpha_mes),2) * std::pow(alpha_mes,2) );
  return value2 - value1;
}

double CharmDecayWidths::ARO2_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1p  = 4.0/(std::pow(alpha_lam,2));
  double value1pp = 1.0/(3.0*std::pow(alpha_mes,2));
  double value1   = 0.25*(value1p + value1pp);
  double value2   = 1.0/(48.0 * std::pow(alpha_mes,4) * std::pow(ARO0_DWAVE(alpha_rho, alpha_mes),2) );
  return std::pow(value1-value2, 0.5);
}

double CharmDecayWidths::BRO1_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  

  double value1  = k_value / ( 4*std::pow(6,0.5) * std::pow(alpha_mes,2) );
  double value2  = k_value / ( 2*std::pow(6,0.5) * alpha_rho * alpha_lam );
  double value3p = k_value / ( 4*std::pow(3,0.5) * std::pow(alpha_mes,2) * std::pow(ARO0_DWAVE(alpha_rho, alpha_mes), 2) );
  double value3pp= 1.0/( 2*std::pow(2,0.5)*alpha_lam*alpha_rho ) + 1.0/(4*std::pow(2,0.5)*alpha_mes*alpha_mes);
  double value3  = value3p*value3pp;
  return value1 + value2 - value3;
}

double CharmDecayWidths::F0TOT_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2 = BRO1_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes) ;
  value2 = std::pow(value2,2);
  double value3 = F0_DWAVE(k_value, alpha_rho, alpha_lam);
  return value1 + value2 + value3; 
}

double CharmDecayWidths::I01B0_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = C0_DWAVE(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_DWAVE(alpha_rho, alpha_mes), 3) );
  double value1pp = BRO1_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 2);
  double value1   = value1p*std::pow(value1pp, 2);
  double value2   = value1pp*ARO1_DWAVE(alpha_rho, alpha_mes) + BRO2_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value3   = std::pow(15,0.5) / (4.0 * std::pow( ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 3) );
  return value1*value2*value3*pi_val*pi_val;
}

double CharmDecayWidths::I02B0_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0_DWAVE(alpha_rho, alpha_lam, alpha_mes) / std::pow(ARO0_DWAVE(alpha_rho, alpha_mes), 3) ;
  double value2 = BRO1_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value) / std::pow(ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 2);
  double value3 = (3*std::pow(10.,0.5) * ARO1_DWAVE(alpha_rho, alpha_mes)) / (8.* std::pow( ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 5));
  return value1*value2*value3*pi_val*pi_val;
}

double CharmDecayWidths::I01B0TOT_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  value2 = -0.0236529;
  return value1*value2;
}

double CharmDecayWidths::I02B0TOT_DWAVE(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}


// RADIAL excitations
// lambda wave functions
double CharmDecayWidths::CBARFIN_RADIAL(double alpha_rho, double alpha_lam){
  double value1 = std::pow(3,0.75) * std::pow(2./3.,0.5);
  double value2 = std::pow(1.0/( std::pow(pi_val, 4) * std::pow(alpha_lam, 2) *std::pow(alpha_rho, 2) ), 0.75);
  return value1 * value2 ;// * mycomplex;//define complex if needed
}

//CBARIN_RADIAL --> same as P-wave (CBARIN)
//CMESON_RADIAL --> same as P-wave (CMESON)

double CharmDecayWidths::C0_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1 = 1.0 / (3.0*std::pow(3,0.5));
  double value2 = CBARFIN_RADIAL(alpha_rho, alpha_lam) * CBARIN(alpha_rho, alpha_lam) * CMESON(alpha_mes);
  return value1 * value2;
}

double CharmDecayWidths::ARO0_RADIAL(double alpha_rho, double alpha_mes){
  return ARO0_DWAVE(alpha_rho, alpha_mes);}

double CharmDecayWidths::F0_RADIAL(double k_value, double alpha_rho, double alpha_lam){
  return F0_DWAVE(k_value, alpha_rho, alpha_lam);
}

double CharmDecayWidths::F01_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1p  = 1.0/(2.0*std::pow(2,0.5)*alpha_lam*alpha_rho) + 1./(4*std::pow(2,0.5)*alpha_mes*alpha_mes);  
  double value1pp = value1p*value1p;  
  double value1ppp= 1./(std::pow(ARO0_RADIAL(alpha_rho, alpha_mes),2));  
  double value1 = value1pp*value1ppp;
  double value2 = 1.0/(8.0*std::pow(alpha_mes,2));
  return std::pow(k_value,2)*(value1-value2);
}

double CharmDecayWidths::BRO2_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = 1.0/(2.0*alpha_lam*alpha_rho) + 1.0/(4.0*alpha_mes*alpha_mes);
  double value2 = 1.0/(std::pow(ARO0_RADIAL(alpha_rho, alpha_mes),2));
  return k_value*(value1*value2 - 2.0);
}

double CharmDecayWidths::ARO1_RADIAL(double alpha_rho, double alpha_mes){
  double value1 = std::pow(2./3., 0.5);
  double value2 = 1./(2.0 * std::pow(6,0.5) * std::pow(ARO0_RADIAL(alpha_rho, alpha_mes),2) * std::pow(alpha_mes,2) );
  return value2 - value1;
}

double CharmDecayWidths::ARO2_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes){
  double value1p  = 4.0/(std::pow(alpha_lam,2));
  double value1pp = 1.0/(3.0*std::pow(alpha_mes,2));
  double value1   = 0.25*(value1p + value1pp);
  double value2   = 1.0/(48.0 * std::pow(alpha_mes,4) * std::pow(ARO0_RADIAL(alpha_rho, alpha_mes),2) );
  return std::pow(value1-value2, 0.5);
}

double CharmDecayWidths::BRO1_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){  
  double value1  = k_value / ( 4*std::pow(6,0.5) * std::pow(alpha_mes,2) );
  double value2  = k_value / ( 2*std::pow(6,0.5) * alpha_rho * alpha_lam );
  double value3p = k_value / ( 4*std::pow(3,0.5) * std::pow(alpha_mes,2) * std::pow(ARO0_RADIAL(alpha_rho, alpha_mes), 2) );
  double value3pp= 1.0/( 2*std::pow(2,0.5)*alpha_lam*alpha_rho ) + 1.0/(4*std::pow(2,0.5)*alpha_mes*alpha_mes);
  double value3  = value3p*value3pp;
  return value1 + value2 - value3;
}

double CharmDecayWidths::F0TOT_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = F01_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value2 = BRO1_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value) / ARO2_RADIAL(alpha_rho, alpha_lam, alpha_mes) ;
  value2 = std::pow(value2,2);
  double value3 = F0_RADIAL(k_value, alpha_rho, alpha_lam);
  return value1 + value2 + value3; 
}

double CharmDecayWidths::I01B0_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0_RADIAL(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_RADIAL(alpha_rho, alpha_mes), 3) );
  double value2p = BRO1_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 2);
  double value2 = value2p*ARO1_RADIAL(alpha_rho, alpha_mes) + BRO2_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value3 = std::pow(12.,0.5)/(4 * std::pow(ARO2_RADIAL(alpha_rho, alpha_lam, alpha_mes) , 3));
  double value4 = 1.5 - std::pow(value2p,2)/std::pow(alpha_lam,2);
  return value1 * value2 * value3 * value4 * std::pow(pi_val, 2.5);
}

double CharmDecayWidths::I02B0_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = C0_RADIAL(alpha_rho, alpha_lam, alpha_mes) / (std::pow(ARO0_RADIAL(alpha_rho, alpha_mes), 3) );
  double value2p = BRO1_DWAVE(alpha_rho, alpha_lam, alpha_mes, k_value)/std::pow(ARO2_DWAVE(alpha_rho, alpha_lam, alpha_mes), 2);
  double value2 = value2p*ARO1_RADIAL(alpha_rho, alpha_mes) + BRO2_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
  double value3 = (3.*std::pow(12.,0.5))/(8 * std::pow(alpha_lam,2) * std::pow(ARO2_RADIAL(alpha_rho, alpha_lam, alpha_mes) , 5) );
  return value1 * value2 * value3 * std::pow(pi_val, 2.5);
}

double CharmDecayWidths::I01B0TOT_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I01B0_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}

double CharmDecayWidths::I02B0TOT_RADIAL(double alpha_rho, double alpha_lam, double alpha_mes, double k_value){
  double value1 = std::exp(F0TOT_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value));
  double value2 = I02B0_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value);
  return value1*value2;
}


#endif

//TEST RADIAL
  // std::cout<<"test mathematica::   CBARFIN  "<<CBARFIN_RADIAL(alpha_rho, alpha_lam)<<std::endl;
  // std::cout<<"test mathematica::   CBARIN  "<<CBARIN(alpha_rho, alpha_lam)<<std::endl;
  // std::cout<<"test mathematica::   CMESON  "<<CMESON(alpha_mes)<<std::endl;  
  // std::cout<<"test mathematica::   C0  "<<C0_RADIAL(alpha_rho, alpha_lam, alpha_mes)<<std::endl;
  // std::cout<<"test mathematica::   ARO0  "<<ARO0_RADIAL(alpha_rho, alpha_mes)<<std::endl;
  // std::cout<<"test mathematica::   F0  "<<F0_RADIAL(k_value, alpha_rho, alpha_lam)<<std::endl;
  // std::cout<<"test mathematica::   F01  "<<F01_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value)<<std::endl;
  // std::cout<<"test mathematica::   BRO2  "<<BRO2_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value)<<std::endl;
  // std::cout<<"test mathematica::   ARO1  "<<ARO1_RADIAL(alpha_rho, alpha_mes)<<std::endl;
  // std::cout<<"test mathematica::   ARO2  "<<ARO2_RADIAL(alpha_rho, alpha_lam, alpha_mes)<<std::endl;
  // std::cout<<"test mathematica::   BRO1  "<<BRO1_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value)<<std::endl;
  // std::cout<<"test mathematica::   F0TOT  "<<F0TOT_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value)<<std::endl;
  // std::cout<<"test mathematica::   I01B0  "<<I01B0_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value)<<std::endl;
  // std::cout<<"test mathematica::   I02B0  "<<I02B0_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value)<<std::endl;
  // std::cout<<"test mathematica::   I01B0TOT  "<<I01B0TOT_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value)<<std::endl;
  // std::cout<<"test mathematica::   I02B0TOT  "<<I02B0TOT_RADIAL(alpha_rho, alpha_lam, alpha_mes, k_value)<<std::endl;
