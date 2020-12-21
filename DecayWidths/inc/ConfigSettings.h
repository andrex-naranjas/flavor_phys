//Config Settings
#ifndef CONFIGSETTINGS_H
#define CONFIGSETTINGS_H

#include <vector>
#include <string>

struct Config{

  std::string InputFileDir;
  std::string OutputFileDir;
  std::string PrintResults;
  std::vector<std::string> ModeExc;
  std::vector<std::string> decayProd;
  std::vector<double> MassA;
  std::vector<double> MassB;
  std::vector<double> MassC;
  std::vector<double> JA_quantum;
  std::vector<double> LA_quantum;
  std::vector<double> SA_quantum;
  std::vector<double> SB_quantum;
  std::vector<double> SC_quantum;
  std::vector<double> Slight_quantum;
  std::vector<double> Slightf_quantum;
  std::vector<double> S_quantum;
  std::vector<double> S1_quantum;
  std::vector<double> S2_quantum;
  std::vector<double> S3_quantum;
  std::vector<double> S4_quantum;
  std::vector<double> S5_quantum;      

};

#endif
