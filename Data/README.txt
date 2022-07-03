The folder Data/PPG_stable contains the concentrations predicted by the Langevin model.
The folder Data/BayOpt contains the improved concentrations that were outputted by the data assimilation system.

The file 'optParamsFixPred.txt' in the Data/Table_5.1 folder contains the optimization results presented in Table 5.1 p.109.
-	Exp. Nr. refers to the corersponding number of the Prairie Grass experiment.
-	wDir refers to the wind direction
-	SDAzim refers to the standard deviation of the wind direction
-	Q refers to the source strength
-	the subindex m refers to the measured value
-	the subindex p refers to the improved value by the data assimilation tool (prediction)
-	fval0 is the value of the cost function for the initial guess
-	fval1 is the value of the cost function after optimization
-	file nr. refers to the number in the file names fixPred<number> and mesoFixPred<number>, which contain the direct output from the data assimilation tool.

The file 'optParamsFixObsFinal.txt' in the Data/Table_5.2 folder contians the optimization results presented in Table 5.2 p. 111.
-	The column names have the same meaning as for file 'optParamsFixPred.txt'
-	file nr. refers to the number in the file names fixObsUnif<number>, mesoFixObs<number> and fixObsNoMeso<number>, which contain the direct output from the data assimilation tool.

The folder Data/inputFiles is an example of how the run directory should be set up in order to run the C++ code.
The files PGARCS.txt and PGVertConc.txt contain the concentration measurements along the arcs and towers in Project Prairie Grass.
The script for running the Langevin model on the Project Prairie Grass data set has been provided in the files 'script_stable.pbs' (stable atmospheric stratification) and
'script_unstable.pbs' (unstable stratification). Here, genInputFiles is the executable from the MATLAB routine modelInputFiles. This routine generates the 
PARAMS_UopenField en PARAMS_KTaylorOF files, which are required as input for the Langevin model. This routine can be found in the folder Code/MATLAB.

The structure of the file 'modelParamList.txt' in the folder Data/inputFiles/case.orig is as follows
<exp. nr.>	<L>		<u_*>	<h_i>	<wind direction>	<\sigma_u>		<\sigma_v>		<\sigma_w>		<\sigma_azim>	<wind speed>
with
-	L the Monin-Obukhov parameter
-	u_* the friction velocity
-	h_i the mixing height
-	\sigma is the standard deviation of the respective wind component
-	\sigma_azim is the standard deviation of the horizontal wind direction
The file 'stationCoord.txt' in case.orig contains the coordinates of all the measurement locations in Project Prairie Grass. 
The content of the file 'meteo5606.txt' in case.orig can be ignored, now it is only used to determine the number of times that the meteodata needs to be refreshed (one in this case).
In the file 'inputDict' in case.orig, only the value of the following inputs should be considered:
-	meteoFile
-	terrainModel
-	stackHeight (meter)
-	releaseDuration (total duration of the pollutant release in seconds)
-	latitude (degree)
-	measureFreq (time in seconds after which the meteo data is refreshed)
-	noPartPerRel (the number of particles that are each time released at once)
-	xMax, yMax, zMax (downwind extension of the domain in meter)
-	sourceStrength (in gram per second)
-	dtByTauL	(ratio of the time step and the Lagrangian time scale)
-	z0	(roughness length in meter)

The data assimilation tool has been implemented in MATLAB. The main file 'bayOpt.m' can be found in the folder Code/MATLAB.
A test case can be run by using the Langevin puff model (LangevinPuff.m).
In the file 'bayOpt.m' the observations variable (line 7) needs then to be set equal to 'predictions'.
The folder PPG_stable in the Data folder is expected to be found in the output folder by the data assimilation tool.
The conc_exp0.txt file in the output.orig directory contains the LangevinPuff concentrations using a wind direction of 180° and a source strength of 100 g/s.
The bayOpt.m can then be run like bayOpt(0,'windDirection',170,'sourceStrength',50) with 0 the experiment number and the other two inputs are the initial guesses.

If a parameter optimization for field measurements needs to be conducted, then observations must equal 'measurements'. 
In the latter case, lines 15, 16, 31, 32, 34, 35 and 41 need to be uncommented and line 42, where LangevinPuff is called, needs to be commented.
For example, the code can then be run as bayOpt(23,'windDirection',128,'windDirStdDev',7.30,'sourceStrength',40.9). 
Here, the numbers 128°, 7.30° and 40.9 g/s represent the measured values that need to be optimized.

The rows in the CONC, JAC and HESS output files correspond with the value at the corresponding station (row) in the file 'stationsCoord.txt'.
The third and fourth column in the JAC file contain the multiplication factors for the concentration kernel estimator in order to calculate the first-order 
partial derivatives of the concentration to the wind direction and its standard deviation respectively.
The third and sixth column in the HESS file contain the multiplication factors for the concentration kernel estimator in order to calculate the second-order 
partial derivatives to the wind direction and its standard deviation respectively. The fourth and fifth column serve for the calculation of the second-order 
mixed partial derivatives. The partial derivatives to the source strength are calculated by the data assimilation tool in MATLAB.
