//Load Settings includes
#ifndef LOADSETTINGS_H
#define LOADSETTINGS_H

#include "ConfigSettings.h"

#include <string>
#include <vector>

class LoadSettings{

 public:
  LoadSettings();
  virtual ~LoadSettings();

  virtual void loadConfig(Config& config, std::string fileName="config.txt");

  //private:

};

#endif
