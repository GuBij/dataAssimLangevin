The README file in the Data folder provides a despcription of how the data assimilation tool should be used.

The C++ code is the implementation of the Langevin particle dispersion model. 
Besides concentrations also the Jacobian and the Hessian, required for the optimization, are outputted.
This code makes use of MPI and the Mersenne Twister (random number generator) from the Intel Math Kernel Library (MKL). 
The code can be compiled as follows: 
1) Go to the folder Code/C++.
2) Run 'make all' (first time) or 'make fastBayesian'.

The routine bayOpt.m in Code/MATLAB is the main file of the data assimilation tool.
The function PGConc in Code/MATLAB has been used to plot the concentration profiles in Fig. 5.1 (p. 114), Fig. 5.2 (p. 115), Fig. 5.3 (p. 117) and Fig. 5.4 (p. 119).
The function predSkill in Code/MATLAB has been used to calculate the Chang and Hanna statistics in Table 5.3 (p. 112).
