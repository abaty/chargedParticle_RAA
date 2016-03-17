#include <string>
#include <iostream>

#include "TString.h"
#include "TObject.h"

class EventSelectionCorrector
{
  public:
    EventSelectionCorrector();
    double getEventWeightFromMC(int M, int nVtx);
    double getEventWeightFromData(int M, int nVtx);
    double getZeroMultFrac();
    virtual ~EventSelectionCorrector();

  private:
  
    static const double trigEff_nVtx1_fromMC[30];
    static const double trigEff_nVtx2_fromMC[30];
    static const double trigEff_nVtx1_fromData[30];
    static const double trigEff_nVtx2_fromData[30];
    static const double zeroMFraction;
};

const double EventSelectionCorrector::trigEff_nVtx1_fromMC[30] = {
        0., 0.88496, 0.982026, 0.993608, 0.997805, 0.998806, 0.998962, 0.9997, 0.999269, 0.99968, 
        0.998194, 1, 0.999696, 0.999657, 0.999188, 0.999548, 1, 1, 1, 0.999234, 
        1, 0.998584, 0.998328, 1, 0.997992, 1, 1, 1, 0.997674, 1 };

const double EventSelectionCorrector::trigEff_nVtx2_fromMC[30] = {
        0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

const double EventSelectionCorrector::trigEff_nVtx1_fromData[30] = {
        0., 0.826584, 0.968855, 0.995067, 0.994862, 0.999544, 0.999634, 0.999151, 0.999752, 0.999713, 
        0.99936, 0.999622, 0.999115, 0.999502, 0.999465, 0.999354, 1, 1, 0.995943, 1, 
        1, 1, 1, 1, 0.997625, 0.997253, 0.996732, 1, 1, 1 };

const double EventSelectionCorrector::trigEff_nVtx2_fromData[30] = {
        0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

const double EventSelectionCorrector::zeroMFraction = 0.;


EventSelectionCorrector::EventSelectionCorrector()
{
}

EventSelectionCorrector::~EventSelectionCorrector()
{
}

double EventSelectionCorrector::getEventWeightFromMC(int M, int nVtx)
{
  if(nVtx==1) {
     if( M==0 || M>29) return 1;
     return 1./trigEff_nVtx1_fromMC[M];
  }
  if(nVtx==2) {
     if( M==0 || M>29) return 1;
     return 1./trigEff_nVtx2_fromMC[M];
  }
  return 1;
}

double EventSelectionCorrector::getEventWeightFromData(int M, int nVtx)
{
  if(nVtx==1) {
     if( M==0 || M>29) return 1;
     return 1./trigEff_nVtx1_fromData[M];
  }
  if(nVtx==2) {
     if( M==0 || M>29) return 1;
     return 1./trigEff_nVtx2_fromData[M];
  }
  return 1;
}

double EventSelectionCorrector::getZeroMultFrac()
{
  std::cerr<<" 'zeroMFraction' NOT SUPPORTED" << std::endl;
  return zeroMFraction;
}

