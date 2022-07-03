%{
delete run_genInputFiles.sh
delete genInputFiles

mcc -v -R -nodisplay -m modelInputFiles.m -o genInputFiles
%}


 delete run_runBayesian.sh 
 delete runBayesian 
 mcc -v -R -nodisplay -m bayOpt.m -o runBayesian 
