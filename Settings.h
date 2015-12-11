#ifndef INPUTSETTINGS
#define INPUTSETTINGS

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

class Settings {
  public:
 
  const int ntrkBins = 42;
  float xtrkbins[ntrkBins+1] = { 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 , 1.1 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 , 3.2 , 4.0 , 4.8 , 5.6 , 6.4 , 7.2 , 9.6 , 12.0, 14.4,19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8,73.6,86.4,103.6,120.8,140,165,190,220,250,280,310,350,400};
  const int nTriggers = 4;

  const int njetBins = 240;
  const int maxJetBin = 1200;

  Settings();
};
 
Settings::Settings()
{
  std::cout << "Getting settings from file: " << inputFile << std::endl;
}

#endif
