#ifndef UhomTurb_H
#define UhomTurb_H

#include "UopenField.h"

 class UhomTurb
 :
   public UopenField
 {

   public:
    UhomTurb(const dictionary&, bool);
    void U ( double* , const double& );
 };

#endif
