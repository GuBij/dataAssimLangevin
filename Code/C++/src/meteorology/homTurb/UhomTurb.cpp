#include "UhomTurb.h"

 UhomTurb::UhomTurb(const dictionary& dict, bool master_proc )
 :
 UopenField(dict, master_proc)
 {}

 void UhomTurb::U(double* u, const double& z)
 {
	UopenField::U(u,UopenField::heff());
 }


