#ifndef DENSITYKERNEL_H
#define DENSITYKERNEL_H

 #include <string>
 using namespace std;

 namespace densityKernel
 {
    void bandWidth
    (
	double**,
	double*,
	int,
	unsigned long,
	bool
    );

    void bandWidth3DKernel
    ( 
	double**, 
	double*, 
	int, 
	unsigned long 
    );

    void bandWidthProdKernel
    (
        double**,
        double*,
        int,
        unsigned long
    );

    void bandWidthProdEpanechnikov
    (
        double**,
        double*,
        int,
        unsigned long
    );

    void bandWidthGaussian
    (
        double**,
        double*,
        int,
        unsigned long
    );

    double weight
    ( 
	string 
    );

    unsigned fac
    (
	unsigned
    );

    double integralKernel
    ( 
	bool 
    );

    double integral3DKernel
    (
    );

   double integralProdKernel
   (
   );

   double integralProdEpanechnikov
   (
   );

   double integral3DEpanechnikov
   (
   );

   double integralGaussian
   (
   );

   void bandWidthHist
   ( 
	double*, 
	double*, 
	int, 
	unsigned long, 
	double*
   );

   void transformParams
   ( 
	double*, 
	double*, 
	int, 
	unsigned long 
   );

   void transformData
   ( 
	double*, 
	double*, 
	double*, 
	double*, 
	unsigned long, 
	int, 
	double*, 
	unsigned long 
   );

   void transformData
   ( 
	double*, 
	double*, 
	double*,
	double*, 
	unsigned long,
	int,
	double*, 
	double*,
	unsigned long 
   );

   double BezierCurve
   ( 
	double, 
	double*, 
	double*, 
	unsigned, 
	double*
   );

   double primitiveBezierCurve
   (
        double,
        double*,
        double*,
        unsigned,
        double*
   );

   void binomialCoeff
   ( 
	double*, 
	unsigned 
   );

   void localRegression
   (
	double*, 
	double, 
	double*, 
	unsigned, 
	unsigned, 
	double
   );

  };

#endif
