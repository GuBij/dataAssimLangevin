 #include <iostream>
 #include <vector>
 #include <algorithm>
 #include <mathimf.h>
 #include "densityKernel.h"

namespace densityKernel
{

 void bandWidth( double** bandWidthVec, double* distrSigmaSlaves, int noProcs, unsigned long dataSizeProc, bool useProdKernel )
 {
    if ( useProdKernel )
	bandWidthProdEpanechnikov(bandWidthVec, distrSigmaSlaves, noProcs, dataSizeProc);
    else
	bandWidth3DKernel(bandWidthVec, distrSigmaSlaves, noProcs, dataSizeProc);
 }

 void bandWidth3DKernel( double** bandWidthVec, double* distrSigmaSlaves, int noProcs, unsigned long dataSizeProc )
 {
    unsigned long nop = 0;
    double varPos[3] = {};
    for (unsigned i = 0; i < noProcs; i++ )
    {
        nop += distrSigmaSlaves[i*dataSizeProc];
        varPos[0] += distrSigmaSlaves[i*dataSizeProc+1]; // x-dir
        varPos[1] += distrSigmaSlaves[i*dataSizeProc+2]; // z-dir
        varPos[2] += distrSigmaSlaves[i*dataSizeProc+3]; // y-dir
    }
    varPos[0] /= nop; // sqrt(varPos[0]/nop);
    varPos[1] /= nop; 
    varPos[2] /= nop; // sqrt(varPos[2]/nop);

    double A = (varPos[0]*varPos[0]*varPos[1]*varPos[1]+varPos[0]*varPos[0]*varPos[2]*varPos[2]+varPos[2]*varPos[2]*varPos[1]*varPos[1])/(varPos[0]*varPos[0]*varPos[1]*varPos[1]*varPos[2]*varPos[2])+0.5*(varPos[0]*varPos[1]+varPos[1]*varPos[2]+varPos[0]*varPos[2])*(varPos[0]*varPos[1]+varPos[1]*varPos[2]+varPos[0]*varPos[2])/(varPos[0]*varPos[1]*varPos[2]*varPos[0]*varPos[1]*varPos[2]);
    varPos[0] = sqrt(varPos[0]);
    varPos[1] = sqrt(varPos[1]);
    varPos[2] = sqrt(varPos[2]); 
    A = 8*varPos[0]*varPos[1]*varPos[2]*2/A;
    A *= 3*15/14*sqrt(acos(-1))*7*7;
    *bandWidthVec[0] = pow(A/nop,1.0/7.0);
    *bandWidthVec[1] = *bandWidthVec[0];
 }
/*
 void bandWidth3DKernel( double** bandWidthVec, double* distrSigmaSlaves, int noProcs, unsigned long dataSizeProc )
 {
    unsigned long nop = 0;
    double varPos = 0;
    for (unsigned i = 0; i < noProcs; i++ )
    {
	nop += distrSigmaSlaves[i*dataSizeProc];
	varPos += distrSigmaSlaves[i*dataSizeProc+1];

	distrSigmaSlaves[i*dataSizeProc] = -1;
	distrSigmaSlaves[i*dataSizeProc+1] = -1;
	distrSigmaSlaves[i*dataSizeProc+2] = -1;
    }
    varPos = sqrt(varPos/nop);

    vector<double> distrSigma(distrSigmaSlaves, distrSigmaSlaves + noProcs*dataSizeProc);
    sort(distrSigma.begin(),distrSigma.end());

    unsigned long m25=floor(0.25*nop), m75=floor(0.75*nop), first = noProcs*dataSizeProc - nop;
    double IQR = (distrSigma[first + m75]-distrSigma[first + m25])/1.34;

    if ( nop > 100 && IQR < varPos && IQR > 0)
      *bandWidthVec[0] = IQR;
    else
      *bandWidthVec[0] = varPos;

    double beta = 0, alpha = 0, sign = -1.0;
    for ( unsigned k = 0; k < 8; k++)
    {
	sign *= -1.0;
	beta += sign/((2*k+1)*fac(k)*fac(7-k));
	if (k < 6)
	  alpha += sign/((2*k+1)*fac(k)*fac(5-k));
    }
    const double PI = acos(-1);
    beta *= 2*PI*fac(6); alpha *= PI*6;
    double A = 32*PI*sqrt(PI)/5*beta/(alpha*alpha);
    *bandWidthVec[0] *= 0.85*pow(A/nop,1.0/7)/100;
//    h *= pow(nop,-1.0/7.0)*weight("bw");

 }
*/
 void bandWidthProdKernel( double** bandWidthVec, double* distrSigmaSlaves, int noProcs, unsigned long dataSizeProc )
 {
    unsigned long nop = 0;
    double varPos[2] = {};
    for (unsigned i = 0; i < noProcs; i++ )
    {
        nop += distrSigmaSlaves[i*dataSizeProc];
        varPos[0] += distrSigmaSlaves[i*dataSizeProc+1];
	varPos[1] += distrSigmaSlaves[i*dataSizeProc+2];
    }
    varPos[0] = sqrt(varPos[0]/nop);
    varPos[1] = sqrt(varPos[1]/nop);

    vector<double> distrSigma(2*nop,0);
    vector<double> distrSigmaVec(distrSigmaSlaves, distrSigmaSlaves + noProcs*dataSizeProc);
    unsigned insertAt = 0;
    for (unsigned i = 0; i < noProcs; i++ )
    {
        unsigned ii = i*dataSizeProc, nopi = long(distrSigmaVec[ii]);
        move(distrSigmaVec.begin()+ii+3,distrSigmaVec.begin()+ii+2+nopi,distrSigma.begin()+insertAt);
        move(distrSigmaVec.begin()+ii+3+(dataSizeProc-3)/2,distrSigmaVec.begin()+ii+2+(dataSizeProc-3)/2+nopi,distrSigma.begin()+nop+insertAt);
        insertAt += distrSigmaVec[ii];
    }
    distrSigmaVec.clear();
    sort(distrSigma.begin(),distrSigma.begin()+nop-1);
    sort(distrSigma.begin()+nop,distrSigma.end());

    unsigned long m25=floor(0.25*nop), m75=floor(0.75*nop);
    double IQR = (distrSigma[m75]-distrSigma[m25])/1.34;
    if ( nop > 100 && IQR < varPos[0] && IQR > 0)
      *bandWidthVec[0] = IQR;
    else
      *bandWidthVec[0] = varPos[0];

    IQR = (distrSigma[nop+m75]-distrSigma[nop+m25]); //1.34;

      double norm2ndDeriv = 0;

    double A = 4*20*20.0/7*4;
    *bandWidthVec[0] *= 0.85*pow(A/nop,1.0/6);

    double beta = 0, alpha = 0, sign = -1.0;
    for ( unsigned k = 0; k < 7; k++)
    {
        sign *= -1.0;
        beta += sign/((2*k+1)*fac(k)*fac(6-k));
        if (k < 4)
          alpha += sign/((2*k+3)*fac(k)*fac(3-k));
    }
    beta *= 2*fac(6); alpha *= 2*6;
    if ( norm2ndDeriv == 0 )
    {
        A = 8*sqrt(acos(-1))/3*beta/(alpha*alpha);
        *bandWidthVec[1] = 0.85*varPos[1]*pow(A/nop,0.2); // 0.85*
    }
    else
    {
        A = beta/(alpha*alpha*norm2ndDeriv);
        *bandWidthVec[1] = pow(A/nop,0.2);
    }
 }

 void bandWidthProdEpanechnikov( double** bandWidthVec, double* distrSigmaSlaves, int noProcs, unsigned long dataSizeProc )
 {
    unsigned long nop = 0;
    double varPos[3] = {};
    for (unsigned i = 0; i < noProcs; i++ )
    {
        nop += distrSigmaSlaves[i*dataSizeProc];
        varPos[0] += distrSigmaSlaves[i*dataSizeProc+1]; // x-dir
        varPos[1] += distrSigmaSlaves[i*dataSizeProc+2]; // z-dir
	varPos[2] += distrSigmaSlaves[i*dataSizeProc+3]; // y-dir
    }
    varPos[0] /= nop; // sqrt(varPos[0]/nop);
    varPos[1] = sqrt(varPos[1]/nop); 
    varPos[2] /= nop; // sqrt(varPos[2]/nop);

    double A = (varPos[0]*varPos[0]+varPos[2]*varPos[2])/(varPos[0]*varPos[0]*varPos[2]*varPos[2])+0.5*(varPos[0]+varPos[2])*(varPos[0]+varPos[2])/(varPos[0]*varPos[2]*varPos[0]*varPos[2]);
    varPos[0] = sqrt(varPos[0]);
    varPos[2] = sqrt(varPos[2]); 
    A = 4*varPos[0]*varPos[2]*2/A; 
    A *= 2*4/3*6*6;
    *bandWidthVec[0] = pow(A/nop,1.0/6);

      *bandWidthVec[1] = varPos[1];

   A = 8*sqrt(acos(-1))*5; // 8*sqrt(acos(-1))/3*3/5*5*5
   *bandWidthVec[1] *= pow(A/nop,0.2);
//   *bandWidthVec[1] = pow(15.0/(1.05e-09*nop),0.2);
 }

 void bandWidthGaussian( double** bandWidthVec, double* distrSigmaSlaves, int noProcs, unsigned long dataSizeProc )
 {
    unsigned long nop = 0;
    double varPos[3] = {};
    for (unsigned i = 0; i < noProcs; i++ )
    {
        nop += distrSigmaSlaves[i*dataSizeProc];
        varPos[0] += distrSigmaSlaves[i*dataSizeProc+1]; // x-dir
        varPos[1] += distrSigmaSlaves[i*dataSizeProc+2]; // z-dir
        varPos[2] += distrSigmaSlaves[i*dataSizeProc+3]; // y-dir
    }
    varPos[0] /= nop; // sqrt(varPos[0]/nop);
    varPos[1] = sqrt(varPos[1]/nop);
    varPos[2] /= nop; // sqrt(varPos[2]/nop);

    vector<double> distrSigma(nop,0);
    vector<double> distrSigmaVec(distrSigmaSlaves, distrSigmaSlaves + noProcs*dataSizeProc);
    unsigned insertAt = 0;
    for (unsigned i = 0; i < noProcs; i++ )
    {
        unsigned ii = i*dataSizeProc, nopi = long(distrSigmaVec[ii]);
        move(distrSigmaVec.begin()+ii+4,distrSigmaVec.begin()+ii+3+nopi,distrSigma.begin()+insertAt);
        insertAt += distrSigmaVec[ii];
    }
    distrSigmaVec.clear();
    sort(distrSigma.begin(),distrSigma.end());

    unsigned long m25=floor(0.25*nop), m75=floor(0.75*nop);
    double IQR = (distrSigma[m75]-distrSigma[m25])/1.34;

    double A = (varPos[0]*varPos[0]+varPos[2]*varPos[2])/(varPos[0]*varPos[0]*varPos[2]*varPos[2])+0.5*(varPos[0]+varPos[2])*(varPos[0]+varPos[2])/(varPos[0]*varPos[2]*varPos[0]*varPos[2]);
    varPos[0] = sqrt(varPos[0]);
    varPos[2] = sqrt(varPos[2]);
    A = 2*varPos[0]*varPos[2]*2/A;
    *bandWidthVec[0] = pow(A/nop,1.0/6);

    if ( nop > 100 && IQR < varPos[1] && IQR > 0)
      *bandWidthVec[1] = IQR;
    else
      *bandWidthVec[1] = varPos[1];

   A = 4.0/3.0; 
   *bandWidthVec[1] *= pow(A/nop,0.2);
 }

 double weight(string type)
 {
   double w;

   if (type == "bw")
   {
     w = 720.0/675675; //beta value
     w *= 3*12006225/(4*acos(-1)*36); //inverse alpha squared value
     w *= 32*pow(acos(-1),1.5)/(3*(2+3));
     w = 0.85*pow(3*w,1.0/7.0);
//     w *= 0.85*pow(nop,-1/7);
   }
   else if (type == "kernel")
   {
     w = 0.25*3*315/(48*acos(-1));
   }

   return w;
 }

 unsigned fac(unsigned n)
 {
   unsigned nfac=1;

   if (n==0)
	return nfac;

   for (unsigned i = 1; i < (n+1); i++)
	nfac *= i;

   return nfac;
 }

 double integralKernel( bool useProdKernel )
 {
   double intK;

   if ( useProdKernel )
     intK = integralProdEpanechnikov();
   else 
     intK = integral3DEpanechnikov();

   return intK;
 }

 double integral3DKernel ( )
 {
   unsigned a = 3; double zLow, zUp;

	zLow = -1;
	zUp = 1;

   double intK = (zUp-zLow)/fac(a+1), zUpPow = zUp, zLowPow = zLow;
   int sign = 1;
   for (unsigned k = 1; k < (a+2); k++)
   {
     sign *= -1;
     intK += sign*(zUpPow-zLowPow)/((2*k+1)*fac(k)*fac(a+1-k));
   }

   intK *= fac(a)*acos(-1);

   return intK;
 }

 double integralProdKernel ( )
 {
   const unsigned a = 3; double zLow, zUp;
   
        zLow = -1;
        zUp = 1;

   double intK = (zUp-zLow)/fac(a), zUpPow = zUp, zLowPow = zLow; 
   int sign = 1; 
   for (int k = 1; k < (a+1); k++)
   {
     sign *= -1;
     intK += sign*(zUpPow-zLowPow)/((2*k+1)*fac(k)*fac(a-k));
   }

   intK *= fac(a)*acos(-1)/(a+1);

   return intK;
 }

 double integralProdEpanechnikov ( )
 {
   return 2.0*acos(-1)/3.0;
 }

 double integral3DEpanechnikov ( )
 {
   return 8.0*acos(-1)/15.0;
 }

 double integralGaussian ( )
 {
   return pow(2*acos(-1),1.5);
 }

 void bandWidthHist( double* h, double* distrSigmaSlaves, int noProcs, unsigned long dataSizeProc, double* minMaxVal)
 {
    unsigned long nop = 0;
    double varPos = 0;
    for (unsigned i = 0; i < noProcs; i++ )
    {
        nop += distrSigmaSlaves[i*dataSizeProc];
        varPos += distrSigmaSlaves[i*dataSizeProc+1];

        distrSigmaSlaves[i*dataSizeProc] = -1;
        distrSigmaSlaves[i*dataSizeProc+1] = -1;
    }

    varPos = sqrt(varPos/nop);
    vector<double> distrSigma(distrSigmaSlaves, distrSigmaSlaves + noProcs*dataSizeProc);
    sort(distrSigma.begin(),distrSigma.end());

    unsigned long m25=floor(0.25*nop), m75=floor(0.75*nop), first = noProcs*dataSizeProc - nop;
    double IQR = distrSigma[first + m75]-distrSigma[first + m25];

    if ( nop > 100 && IQR > 0)
    {
      double cube = pow(nop,-1.0/3);
      h[0] = 3.49*varPos*cube; h[1] = h[0];
      h[2] = 2*IQR*cube;
    }
    else
    {
	unsigned noBins = ceil(2*pow(nop,1.0/3));
	h[0] = (minMaxVal[1]-minMaxVal[0])/noBins;
	h[1] = (minMaxVal[3]-minMaxVal[2])/noBins;
	h[2] = minMaxVal[5]/noBins;
    }
 }

 void transformParams( double* params, double* dataContainer, int noRecords, unsigned long dataSize ) //double* range
 {
    // assume noRecords = 1; code should be modified for noRecords > 1
    vector<double> sortedData(dataContainer, dataContainer + noRecords*dataSize);
    sort(sortedData.begin(),sortedData.end());

    unsigned long nop = dataSize;
    unsigned long m25=floor(0.25*nop), m75=floor(0.75*nop);
    double IQR = sortedData[m75]-sortedData[m25];
    double A = 8*sqrt(acos(-1))*5;

    double binSize = 2*IQR*pow(nop,-1.0/3.0);
    double upperBnd = ceil(sortedData.back());
    const unsigned noBins = ceil(upperBnd/binSize); binSize = upperBnd/noBins;

    double* pdfApprox = new double[noBins]();
    for ( unsigned long i = 0; i < nop; i++)
    {
      unsigned ii = floor(sortedData[i]/binSize);
      pdfApprox[ii] += 1.0;
    }

    double b0 = (1.5*pdfApprox[0]-0.5*pdfApprox[1])/(nop*binSize);
    double b1 = (pdfApprox[1]-pdfApprox[0])/(nop*binSize*binSize);
    delete[] pdfApprox;

    if ( b0 <= 0 && b1 < 0 )
    {
	b0 = 0; b1 = 0;
    }
    else if ( b0 < 0 )
	b0 = 0;

    cout << "\n b0: " << b0 << ", b1: " << b1 << ", binSize: " << binSize << endl;

    unsigned long m50=floor(0.5*nop);
    params[2] = 1.0/sortedData[m50];
    params[1] = -b1*sortedData[m50];
    params[0] = b0 - params[1];
 }

 void transformData( double* data, double* jacobian, double* meanDrift, double* grid, unsigned long dataSize, int noRecords, double* params, unsigned long gridSize )
 {
   unsigned r;
   for ( unsigned long i = 0; i < noRecords*dataSize; i++ )
   {
     unsigned long ii = i % dataSize;
     if ( ii == 0 )
     {
        r = i/dataSize; meanDrift[4*r+2] = 0;
     }
     data[i] = params[0]*data[i] - params[1]/params[2]*(exp(-params[2]*data[i])-1); 
     meanDrift[4*r+2] += data[i];
   }

   for ( unsigned long i = 0; i < noRecords*gridSize; i++ )
   {
     unsigned long ii = i % gridSize;
     if ( ii == 0 )
        r = i/gridSize;

     jacobian[i] = params[0] + params[1]*exp(-params[2]*grid[3*i+2]); 
     grid[3*i+2] = params[0]*grid[3*i+2] - params[1]/params[2]*(exp(-params[2]*grid[3*i+2])-1);
   }
 }

 void transformData( double* data, double* jacobian, double* meanDrift, double* grid, unsigned long dataSize, int noRecords, double* controlpnts, double* range, unsigned long gridSize )
 {
   unsigned r;
   double normFactor;
   double* nplusplusC = NULL;
   for ( unsigned long i = 0; i < noRecords*dataSize; i++ )
   {
     unsigned long ii = i % dataSize; 
     if ( ii == 0 )
     {
	r = i/dataSize; meanDrift[4*r+2] = 0; 
	unsigned n = range[3*r+2];
	nplusplusC = new double[n+2](); binomialCoeff(nplusplusC,n+1);
	normFactor = primitiveBezierCurve( range[3*r+1], nplusplusC, controlpnts, range[3*r+2], range+3*r );
     }
     data[i] = primitiveBezierCurve( data[i], nplusplusC, controlpnts, range[3*r+2], range+3*r );
     data[i] /= normFactor;
     if ( data[i] < 0 )
     {
	cout << "\n data[" << i << "]: " << data[i] << endl;
     }
     meanDrift[4*r+2] += data[i];
   }

   double* nC = NULL;   
   for ( unsigned long i = 0; i < noRecords*gridSize; i++ )
   {
     unsigned long ii = i % gridSize;
     if ( ii == 0 )
     {
	r = i/gridSize;
	unsigned n = range[3*r+2];
	nC = new double[n+1](); binomialCoeff(nC,n);
	delete[] nplusplusC;
	nplusplusC = new double[n+2](); binomialCoeff(nplusplusC,n+1);
	normFactor = primitiveBezierCurve( range[3*r+1], nplusplusC, controlpnts, range[3*r+2], range+3*r );
     }
     jacobian[i] = BezierCurve( grid[3*i+2], nC, controlpnts, range[3*r+2], range+3*r );
     grid[3*i+2] = primitiveBezierCurve( grid[3*i+2], nplusplusC, controlpnts, range[3*r+2], range+3*r );
     grid[3*i+2] /= normFactor;
   }
   delete[] nC; delete[] nplusplusC;
 }

 double BezierCurve( double t, double* binomCoeff, double* controlpnts, unsigned n, double* range )
 {
   double* monos = new double[n+1];
   monos[0] = 1;
   for ( unsigned j = 1; j <= n; j++ )
   {
     monos[j] = monos[j-1]*(t-range[0])/(range[1]-range[0]);
   }
   double monos2 = 1, B = 0;
   for ( unsigned k = n; k > 0; k-- )
   {
     B += controlpnts[k]*monos[k]*monos2*binomCoeff[k];
     monos2 *= (range[1]-t)/(range[1]-range[0]);
   }
   delete[] monos;

   return B;
 }

 double primitiveBezierCurve( double t, double* binomCoeff, double* controlpnts, unsigned n, double* range )
 {
   n++;
   double* monos = new double[n+1];
   monos[0] = 1;
   for ( unsigned j = 1; j <= n; j++ )
   {
     monos[j] = monos[j-1]*(t-range[0])/(range[1]-range[0]);
   }
   double primitive = 0;
   double monos2 = 1, B = 0;
   for ( unsigned k = n; k > 1; k-- )
   {
     primitive += monos[k]*monos2*binomCoeff[k]/n;
     B += controlpnts[k-1]*primitive*(range[1]-range[0]);
     monos2 *= (range[1]-t)/(range[1]-range[0]);
   }
   delete[] monos;

   return B;
 }

 void binomialCoeff( double* coeffVec, unsigned n )
 {
   coeffVec[0] = 1;
   unsigned ii = 1, k = 2;
   double save[2] = {1,0};
   for ( unsigned long i = 1; i < ((n+1)*(n+2)/2-n); i++ )
   {
     if ( ii == k )
     {
	coeffVec[ii-1] = save[0]; ii = 1; k++; 
     }
     save[1] = coeffVec[ii-1]+coeffVec[ii];
     coeffVec[ii-1] = save[0];
     save[0] = save[1];
     ii++;
   }
   coeffVec[n] = save[0];
 }

 void localRegression(double* y, double meshWidth, double* samples, unsigned noGridPoints, unsigned noSamples, double h)
 {
   unsigned startIndex = 0;
   for (unsigned i = 0; i < noGridPoints; i++)
   {
     const double xi = (i+0.5)*meshWidth;
     int k = startIndex-1;
     bool outside = samples[startIndex] <= (xi-h);
     double xlocalbar, ylocalbar, sumWeight;
     while ( ( k+1 <= noSamples) && (outside || abs(xi - samples[k+1]) < h ))
     {
       k = k + 1;
       if ( outside )
        startIndex = k;
       else
       {
	 double weight = 1-(xi-samples[k])*(xi-samples[k])/(h*h);
	 weight = weight*weight*weight;
	 xlocalbar += weight*samples[k];
	 ylocalbar += weight*samples[k+noSamples];
	 sumWeight += weight;
       }
       outside = samples[k,1] <= (xi-h);
     }
     xlocalbar /= sumWeight;
     ylocalbar /= sumWeight;
     double cov = 0, var_x = 0;
     for (unsigned j = startIndex; j < (k+1); j++)
     {
	cov += (samples[j]-xlocalbar)*(samples[j+noSamples]-ylocalbar);
	var_x += (samples[j]-xlocalbar)*(samples[j]-xlocalbar);
     }
     y[i] = ylocalbar + cov/var_x*(xi-xlocalbar);
   }
 }

}
