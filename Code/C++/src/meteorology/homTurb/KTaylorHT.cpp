#include "KTaylorHT.h"

 KTaylorHT::KTaylorHT(const dictionary& dict, bool master_proc )
 :
 KTaylorOF(dict, master_proc)
 {}

 void KTaylorHT::eddyDiff ( double* eddyVec, const double& z)
 {
    double tL[3];
    KTaylorOF::tauL(tL,0);

    *eddyVec=sqrt(2/tL[0])*sigma_u;
    *(eddyVec+1)=sqrt(2/tL[1])*sigma_v;
    *(eddyVec+2)=sqrt(2/tL[2])*sigma_w;
 }

 double KTaylorHT::nonHomTurbCorr ( const double& z)
 {
	return 0.0;
 }

