#ifndef KTaylorHT_H
#define	KTaylorHT_H

 #include "KTaylorOF.h"

 class KTaylorHT
 :
 public KTaylorOF
 {

  public:
   KTaylorHT(const dictionary&, bool);
   void eddyDiff( double*, const double& );
   double nonHomTurbCorr( const double& );
 };

#endif
