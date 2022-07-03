#include "KTaylorOF.h"

 KTaylorOF::KTaylorOF(const dictionary& dict, bool master_proc )
 :
 openField(dict,true, master_proc),
 name_("KTaylorOF"),
 dtByTauL_(dict.getfl(11)),
// minTimeStep_(dict.getfl(12)),
 mfreq(dict.getfl(7))
 {
    eddyParams = new double[noRecords*8];
    ifstream eddyParamFile("output/PARAMS_" + name_ + "_" + dict.get(0));
    for (unsigned i = 0; i < noRecords; i++)
    {
        unsigned timeIndex, ii = 8*i;
        eddyParamFile >> timeIndex >> eddyParams[ii] >> eddyParams[ii+1] >> eddyParams[ii+2] >> eddyParams[ii+3] >> eddyParams[ii+4] >> eddyParams[ii+5] >> eddyParams[ii+6] >> eddyParams[ii+7];
    }
    eddyParamFile.close();

    w2bar = new double[noRecords*noGridPointsZ]();
    w2barGrad = new double[noRecords*noGridPointsZ]();
    tauL_w = new double[noRecords*noGridPointsZ]();
    tauL_u = new double[noRecords*noGridPointsZ]();
    corr = new double[2*noRecords*noGridPointsZ]();
    varianceRate = new double[2*noRecords*noGridPointsZ]();
    timeStep_ = new double[noRecords*noGridPointsZ]();
    timeStep_u = new double[noRecords*noGridPointsZ]();
    double pct, lengthScale, tauL_uv = 0.0;
    bool neutral, updateTimeStep, update = true;
    unsigned heffCell;
    unsigned long i = 0;
    while ( i < (noRecords*noGridPointsZ) )
    {
        unsigned long ii = i % noGridPointsZ, r = i/noGridPointsZ;

        if ( update )
        {
          updateParams(r); unstable = (L < 0 && L > -500); neutral = (abs(L) >= 500); stable = (L>0 && L<500); updateTimeStep = true; update = false; pct = dtByTauL_; heffCell = heff_/resZ;
	  if ( heffCell > (noGridPointsZ-1) )
	     heffCell = noGridPointsZ-1;

	  i += heffCell; ii = heffCell;

	  if ( unstable )
	  {
	    tauL_uv = 0.15*hABL/sigma_u;
	    minTimeStep_ = 10*pct*0.1*30*z0/(propConst*pow(30*z0,1.0/3.0)*(0.55+0.38*29*z0/L));
	  }
	  else if ( stable )
	    minTimeStep_ = pct*0.1*hABL*pow(30*z0/hABL,0.8)/sigma_w;
	  else
	    minTimeStep_ = pct*0.5*30*z0/sigma_w/(1+15*fCori*30*z0/uStar);

        }
	double zi = resZ*(ii+0.5);

	if ( neutral )
	{
	  w2bar[i] = sigma_w;
          tauL_w[i] = 0.5*zi/sigma_w/(1+15*fCori*zi/uStar);
	}
	else if ( unstable )
	{
	  if ( zi < 30*z0 )
	  {
	   w2bar[i] = propConst*pow(30*z0,1.0/3.0);
	   w2barGrad[i] = 0;
	  }
	  else
	  {
	   w2bar[i] = propConst*pow(zi,1.0/3.0);
	   w2barGrad[i] = propConst*propConst*2.0/3.0*pow(1/zi,1.0/3.0);
	  }

	  if ( zi < 0.1*hABL && (zi-z0) < -L )
	  {
	     tauL_w[i] = 0.1*zi/(w2bar[i]*(0.55+0.38*(zi-z0)/L)); 
	  }
	  else if ( zi < 0.1*hABL && (zi-z0) >= -L )
	     tauL_w[i] = 0.59*zi/w2bar[i];
	  else
	     tauL_w[i] = 0.15*hABL/w2bar[i]*(1-exp(-5*zi/hABL));

          if ( tauL_w[i] > tauL_uv ) // ensure proper functioning code; normally occurs when z > hABL, model not valid there
            tauL_w[i] = tauL_uv;
	}
	else 
	{
	  w2bar[i] = sigma_w;
	  tauL_u[i] = 0.085*sqrt(hABL*zi)/sigma_v; //use sigma_u = 1.5*sigma_v and assume tauL_u = tauL_v
	  tauL_w[i] = 0.1*hABL*pow(zi/hABL,0.8)/sigma_w;

	  if ( tauL_w[i] > tauL_u[i] ) // ensure proper functioning code; normally occurs when z > hABL, model not valid there
		tauL_w[i] = tauL_u[i];
	}

	double timeStepi, timeStepi_u; 
	double tauL_wi, tauL_ui;
	if ( stable )
	{
	   tauL_wi = tauL_w[i];
	   tauL_ui = tauL_u[i];
	   timeStepi = pct*tauL_wi;
	   timeStepi_u = pct*tauL_ui;
	}
	else // if ( unstable )
	{
	   tauL_wi = tauL_w[i];
	   timeStepi = pct*tauL_wi;
	   if ( unstable )
		timeStepi *= 10.0;
//	   if ( pct > pct_MIN && pct/w2barGrad[i] < timeStepi )
//		timeStepi = pct/w2barGrad[i];
	}

	if ( stable && timeStepi_u < minTimeStep_ )
	{
           timeStep_u[i] = minTimeStep_;
	}
	else if ( stable )
	{
	   timeStep_u[i] = timeStepi_u;
	}

	if ( timeStepi < minTimeStep_ )
	{
           timeStep_[i] = minTimeStep_;
	}
	else
	   timeStep_[i] = timeStepi;

	if ( ii == heffCell )
        {

	  timeStep_[i] = mfreq/ceil(mfreq/timeStep_[i]); timeStepi = timeStep_[i];
	  if ( timeStep_[i] > minTimeStep_ )
	    pct = timeStep_[i]/tauL_wi;

	  if ( stable &&  timeStep_[i] > minTimeStep_ )
	  {
	    timeStep_u[i] = mfreq/ceil(mfreq/timeStep_u[i]); pct = timeStep_u[i]/tauL_ui; timeStepi_u = timeStep_u[i];
	  }
	  else if ( stable )
	  {
	    timeStep_u[i] = mfreq/ceil(mfreq/timeStep_u[i]); timeStepi_u = timeStep_u[i];
	  }
        }

	if ( neutral )
	{
	  corr[i] = exp(-timeStep_[i]/tauL_wi);
	  corr[i+noRecords*noGridPointsZ] = corr[i];
	  varianceRate[i] = sqrt(2/tauL_wi);
	  varianceRate[i+noRecords*noGridPointsZ] = varianceRate[i]*w2bar[i];
	}
	else if ( unstable )
	{
	  double dt = 0.15*hABL/sigma_u*dtByTauL_; dt = mfreq/ceil(mfreq/dt);
	  corr[i] = exp(-dt*sigma_u/(0.15*hABL));
	  corr[i+noRecords*noGridPointsZ] = exp(-timeStep_[i]/tauL_wi);
	  varianceRate[i] = sqrt(2*sigma_u/(0.15*hABL));
	  varianceRate[i+noRecords*noGridPointsZ] = sqrt(2/tauL_wi)*w2bar[i];
	}
	else
	{
	  corr[i] = exp(-timeStep_u[i]/tauL_ui);
	  corr[i+noRecords*noGridPointsZ] = exp(-timeStep_[i]/tauL_wi);
	  varianceRate[i] = sqrt(2/tauL_ui);
	  varianceRate[i+noRecords*noGridPointsZ] = sqrt(2/tauL_wi)*w2bar[i];
	}

	i++;

	if ( ii == heffCell && heffCell == 0 )
        {
	  updateTimeStep = false;
	}
	else if ( ii == heffCell )
        {
	  i = r*noGridPointsZ;
	}
	else if ( heffCell > 0 && ii == heffCell - 1 )
        {
          i++; updateTimeStep = false;
        }

	if ( i == (r+1)*noGridPointsZ )
	{
	  update = true;
	}
    } 

    updateParams(0);
 }

 void KTaylorOF::tauL ( double* tauVec, const double& z )
 {
        unsigned cellIndex = z/resZ;
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

        cellIndex = cellIndex + metId*noGridPointsZ;
	if ( stable )
	{
	   *tauVec = tauL_u[cellIndex];
	   *(tauVec+1) = *tauVec;
	   *(tauVec+2) = tauL_w[cellIndex];
	}
        else if ( unstable )
        {
           *tauVec = 0.15*hABL/sigma_u;
           *(tauVec+1) = *tauVec;
           *(tauVec+2) = tauL_w[cellIndex];
        }
	else
	{
           *tauVec = tauL_w[cellIndex];
           *(tauVec+1) = tauL_w[cellIndex];
           *(tauVec+2) = tauL_w[cellIndex];
	}
 }

 void KTaylorOF::eddyDiff ( double* eddyVec, const double& z)
 {
        unsigned cellIndex = z/resZ;
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

	cellIndex = cellIndex + metId*noGridPointsZ;

        *eddyVec=varianceRate[cellIndex]*sigma_u;
        *(eddyVec+1)=varianceRate[cellIndex]*sigma_v;
        *(eddyVec+2)=varianceRate[cellIndex+noRecords*noGridPointsZ];
 }

 void KTaylorOF::eddyDiffRough ( double* eddyVec, const double* tauLVec, const double& z)
 {
        unsigned cellIndex = z/resZ;
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

        cellIndex += metId*noGridPointsZ;

        *eddyVec=sqrt(2/tauLVec[0]*(sigma_u*sigma_u-uStar*uStar*uStar*uStar/(w2bar[cellIndex]*w2bar[cellIndex])));
        *(eddyVec+1)=sqrt(2/tauLVec[1])*sigma_v;
        *(eddyVec+2)=sqrt(2/tauLVec[2])*uStar*uStar/w2bar[cellIndex];
 }

 void KTaylorOF::stdevWind ( double* stdevVec, const double& z)
 {
        unsigned cellIndex = z/resZ;
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

        cellIndex = cellIndex + metId*noGridPointsZ;

	*stdevVec=sigma_u;
	*(stdevVec+1)=sigma_v;
	*(stdevVec+2)=w2bar[cellIndex];
 }

 void KTaylorOF::nonHomTurbCorr ( double& a, const double& z, const double& zdot)
 {
        unsigned cellIndex = z/resZ;
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

        cellIndex = cellIndex + metId*noGridPointsZ;

	if ( unstable )
	  a = 0.5*(zdot*zdot/(w2bar[cellIndex]*w2bar[cellIndex])+1)*w2barGrad[cellIndex];
	else
	  a = 0;
 }

 void KTaylorOF::timeStep( double* timeStepVec, const double& z )
 {
        unsigned cellIndex = z/resZ;
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

	cellIndex = cellIndex + metId*noGridPointsZ;

	*(timeStepVec+1) = timeStep_[cellIndex];
	if ( unstable && timeStep_[cellIndex] < 0.15*hABL/sigma_u*dtByTauL_ )
	{
	   *timeStepVec = 0.15*hABL/sigma_u*dtByTauL_; timeStepVec[0] = mfreq/ceil(mfreq/timeStepVec[0]);
	}
	else if ( stable )
	{
	   *timeStepVec = timeStep_u[cellIndex];
	}
	else
	   *timeStepVec = timeStep_[cellIndex];
 }

 void KTaylorOF::autoCorr( double* corrVec, const double& z )
 {
        unsigned cellIndex = z/resZ;
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

        cellIndex = cellIndex + metId*noGridPointsZ;

	corrVec[0] = corr[cellIndex];
	corrVec[1] = corr[cellIndex+noRecords*noGridPointsZ];
 }

 void KTaylorOF::autoCorr( double* corrVec, const double& z, const double& timeStepVert ) //doesn't work yet for stable
 {
        unsigned cellIndex = z/resZ;
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

        cellIndex = cellIndex + metId*noGridPointsZ;

        corrVec[0] = corr[cellIndex];
        corrVec[1] = exp(-timeStepVert/tauL_w[cellIndex]);
 }

 void KTaylorOF::updateParams ( unsigned meteo_id )
 {
        if (meteo_id < noRecords)
        {
           sigma_u = eddyParams[meteo_id*8];
           sigma_v = eddyParams[meteo_id*8+1];
           sigma_w = eddyParams[meteo_id*8+2];
	   propConst = eddyParams[meteo_id*8+3];
           uStar = eddyParams[meteo_id*8+4];
           L = eddyParams[meteo_id*8+5];
	   heff_ = eddyParams[meteo_id*8+6];
	   hABL = eddyParams[meteo_id*8+7];
	   metId = meteo_id;
	   unstable = ( L < 0 && L > -500 ); 
	   stable = ( L > 0 && L < 500 );
       }
        else
        {
           cout << "\n\t ERROR: meteo record with id " << meteo_id << " does not exist.\n" << endl; 
           exit(EXIT_FAILURE);
        }	
 }
