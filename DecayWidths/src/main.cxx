#include "CharmDecayWidths.h"

#include "ConfigSettings.h"
#include "LoadSettings.h"

#include <string>
#include <vector>
#include <iostream>


int main(int argc, char** argv){

  if(argc!=2){
    std::cout<<"Please provide a config file, you provided "<<argc-1<<" files/arguments"<<std::endl;
    std::cout<<"Try again. BYE!"<<std::endl;
    return 1;
  }


  //Get the cofiguration from configFile.txt
  LoadSettings *settings = new LoadSettings();
  Config config;
  settings->loadConfig(config, argv[1]);

  CharmDecayWidths *m_decays = new CharmDecayWidths();
  m_decays ->execute(config);

  
  return 0;
}
