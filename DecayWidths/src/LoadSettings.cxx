//Simple class to parse arguments from config file
#ifndef LOADSETTINGS_CXX
#define LOADSETTINGS_CXX

#include "LoadSettings.h"
#include "ConfigSettings.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <string>


LoadSettings::LoadSettings()
{
}

LoadSettings::~LoadSettings(){}

void LoadSettings::loadConfig(Config& config, std::string fileName){

  std::ifstream fin(fileName);
  std::string line;
  while (getline(fin,line)) {
    std::istringstream sin(line.substr(line.find("=") + 1));    
    if (line.find("InputFileDir") != -1)
      sin >> config.InputFileDir;
    else if (line.find("OutputFileDir") != -1)
      sin >> config.OutputFileDir;
    else if (line.find("PrintResults") != -1)
      sin >> config.PrintResults;
    else if (line.find("ModeExc") != -1){
      std::string name="";
      while(sin >> name)
	config.ModeExc.push_back(name);
    }else if (line.find("DecayProd") != -1){
      std::string name="";
      while(sin >> name)
	config.decayProd.push_back(name);
    }else if (line.find("MassA") != -1){
      double number=0.;
      while(sin >> number)
	config.MassA.push_back(number);
    }else if (line.find("MassB") != -1){
      double number=0.;
      while(sin >> number)
	config.MassB.push_back(number);
    }else if (line.find("MassC") != -1){
      double number=0.;
      while(sin >> number)
	config.MassC.push_back(number);
    }else if (line.find("JA_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.JA_quantum.push_back(number);
    }else if (line.find("LA_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.LA_quantum.push_back(number);
    }else if (line.find("SA_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.SA_quantum.push_back(number);
    }else if (line.find("SB_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.SB_quantum.push_back(number);
    }else if (line.find("SC_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.SC_quantum.push_back(number);
    }else if (line.find("Slight_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.Slight_quantum.push_back(number);
    }else if (line.find("Slightf_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.Slightf_quantum.push_back(number);
    }else if (line.find("S_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.S_quantum.push_back(number);
    }else if (line.find("S1_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.S1_quantum.push_back(number);
    }else if (line.find("S2_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.S2_quantum.push_back(number);
    }else if (line.find("S3_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.S3_quantum.push_back(number);
    }else if (line.find("S4_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.S4_quantum.push_back(number);
    }else if (line.find("S5_quantum") != -1){
      double number=0.;
      while(sin >> number)
	config.S5_quantum.push_back(number);
    }    
    
  }//firstwhile
}//loadConfig

#endif

    // else if (line.find("DoOnlyPlots") != -1)
    //    sin >> config.DoOnlyPlots;
