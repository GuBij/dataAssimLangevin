function [fval,grad,hess]=objFuncBay(x)

 global bckGrndCovMatInv priorWeights y_o x_o sigmaAzimMin varPresent concFile sensFile hessFile inputFileU inputFileK inputFile rmStations modelParams inputFileUOrig inputFileKOrig 

 varID = cumsum(varPresent);

 jacNoCol = 4;
 hessNoCol = 6;
 noControlVar = length(x);
 prior = zeros(noControlVar,1);
 gradPrior = zeros(noControlVar,1);
 HessPrior = zeros(noControlVar);
 keepHessCol = zeros(1,4);
 if ( varPresent(1) )
%  copyfile(sprintf('meteoFiles/%s',inputFileUOrig),sprintf('output/%s',inputFileU));
%  system(strcat('sed -i "s/<wDir>/',num2str(x(1)),'/" output/',inputFileU));
  prior(1) = priorWeights(1)*(x(1)-x_o(1))^2;
  gradPrior(1) = priorWeights(1)*2*(x(1)-x_o(1));
  HessPrior(1,1) = priorWeights(1)*2;
  keepHessCol(1,(1:3)) = 1;
 end
 if ( varPresent(2) ) 
   prior(varID(2)) = priorWeights(varID(2))*(x(varID(2))-sigmaAzimMin); 
   gradPrior(varID(2)) = priorWeights(varID(2));
   keepHessCol(1,(2:4)) = 1;
   sigmaAzimMeso = mesoScaleVar(x(varID(2)),modelParams); 
 else
   sigmaAzimMeso = 0;
   keepHessCol(1,(2:4)) = 0;
 end
% system(strcat('sed -i "s/<mesoVar>/',num2str(sigmaAzimMeso),'/" output/',inputFileU));
% system(sprintf('cp meteoFiles/%s ./output/%s',inputFileKOrig,inputFileK));
 if ( varPresent(3) )
%  system(sprintf('cp case.orig/%s ./',inputFile));
%  system(strcat('sed -i "s/<Q>/',num2str(x(varID(3))),'/" ./',inputFile)); 
  prior(varID(3)) = priorWeights(varID(3))*x(varID(3)); 
  gradPrior(varID(3)) = priorWeights(varID(3));
  keepHessCol = [keepHessCol(1:2),keepHessCol(1),keepHessCol(3:4),keepHessCol(4),keepHessCol(1),keepHessCol(4),1];
 end
 keepHessCol = logical(keepHessCol);
% system('mpirun -np 1 fastPM/fastBayesian -Langevin'); 
 LangevinPuff(x(1),x(varID(3))); % /1000);

 y_p = readFile(concFile,0,3); y_p = y_p(:,3); %predictions
 J_p = readFile(sensFile,0,jacNoCol); J_p = J_p(:,3:end);
 H_p = readFile(hessFile,0,hessNoCol); H_p = H_p(:,3:end);

 y_p(rmStations) = []; J_p(rmStations,:) = []; H_p(rmStations,:)=[];
 if ( sigmaAzimMeso > 0 )
  J_p(:,2) = J_p(:,2)*x(varID(2))/sigmaAzimMeso;
  H_p(:,[2,3]) = H_p(:,[2,3])*x(varID(2))/sigmaAzimMeso;
  H_p(:,4) = H_p(:,4)*(x(varID(2))/sigmaAzimMeso)^2 + J_p(:,2)/x(varID(2));
 end
 if ( varPresent(3) )
  ii = varID(3);
  J_p = [J_p,y_p/x(ii)];
  H_p = [H_p(:,1:2),J_p(:,1)/x(ii),H_p(:,3:4),J_p(:,2)/x(ii),J_p(:,1:2)/x(ii),zeros(size(J_p,1),1)];
 end
 J_p = J_p(:,varPresent);
 H_p = H_p(:,keepHessCol);

 fval = (y_p-y_o)'*bckGrndCovMatInv*(y_p-y_o) + sum(prior); 
 grad = 2*J_p'*bckGrndCovMatInv*(y_p-y_o) + gradPrior;

 hess = zeros(noControlVar);
 for i = 1:noControlVar
  hess(i,:) = 2*J_p(:,i)'*bckGrndCovMatInv*J_p + 2*(y_p-y_o)'*bckGrndCovMatInv*H_p(:,(i-1)*noControlVar+(1:noControlVar));
 end
 hess = hess + HessPrior;

end
