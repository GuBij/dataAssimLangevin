function hess=hessFunc(x,lambda)

 global bckGrndCovMatInv priorCovMatInv y_o x_o concFile sensFile hessFile rmStations modelParams lb ub Qid % priorCovMatInv  

 x = lb + x.*(ub-lb);
 xx_o = lb + x_o.*(ub-lb);

 [sigmaAzimMeso,uStar,U] = mesoScaleVar(x(2),modelParams,x(3));
 noControlVar = length(x);
 minOne=0;
 if ( noControlVar == Qid )
  minOne=1;
 end
 y_p = readFile(concFile,0,3); y_p = y_p(:,3); %predictions
 J_p = readFile(sensFile,0,2+noControlVar-minOne); J_p = J_p(:,3:(2+noControlVar-minOne)); %Jacobian matrix of the predictions
 H_p = readFile(hessFile,0,2+(noControlVar-minOne)^2); H_p = H_p(:,3:(2+(noControlVar-minOne)^2)); %2nd order sensitivities of the predictions

 y_p(rmStations) = []; J_p(rmStations,:) = []; H_p(rmStations,:) = [];
 J_p(:,3) = J_p(:,3)*uStar;
 H_p(:,9) = H_p(:,9)*uStar^2;
 H_p(:,3) = H_p(:,3)*uStar;
 if ( sigmaAzimMeso > 0 ) 
  DiffMesoToTurbVar = x(3)*uStar^2/sigmaAzimMeso*(180/(U*pi))^2;
  J_p(:,3) = J_p(:,3) + J_p(:,2)*DiffMesoToTurbVar;
  H_p(:,9) = H_p(:,9) + H_p(:,8)*DiffMesoToTurbVar*uStar + J_p(:,2)*(1-x(3)*uStar*(180/(sigmaAzimMeso*U*pi))^2)*(180/(U*pi))^2/sigmaAzimMeso*uStar^2;
  H_p(:,3) = H_p(:,3) + H_p(:,2)*DiffMesoToTurbVar; 
  H_p(:,6) = H_p(:,6)*uStar*x(2)/sigmaAzimMeso + H_p(:,5)*DiffMesoToTurbVar*x(2)/sigmaAzimMeso - J_p(:,2)*DiffMesoToTurbVar*x(2)/sigmaAzimMeso^2; 
  J_p(:,2) = J_p(:,2)*x(2)/sigmaAzimMeso;
  H_p(:,[2,4]) = H_p(:,[2,4])*x(2)/sigmaAzimMeso;
  H_p(:,5) = H_p(:,5)*(1-(x(2)/sigmaAzimMeso)^2)/sigmaAzimMeso;
 end
 H_p(:,7) = H_p(:,3);
 H_p(:,8) = H_p(:,6);

 %source strength optimization
 HessPrior = [[2*priorCovMatInv,zeros(1,noControlVar-1-minOne)];zeros(noControlVar-1-minOne,noControlVar-minOne)];
 if ( noControlVar == Qid )
  J_p = [J_p,y_p/x(Qid)];
  H_p = [H_p(:,1:(Qid-1)),J_p(:,1)/x(Qid),H_p(:,Qid:(2*(Qid-1))),J_p(:,2)/x(Qid),H_p(:,(2*(Qid-1)+1):(3*(Qid-1))),J_p(:,3)/x(Qid),J_p(:,1:(Qid-1))/x(Qid),zeros(size(J_p,1),1)];
  HessPrior = [[2*priorCovMatInv,zeros(1,noControlVar-1)];zeros(noControlVar-1,noControlVar)];
 end

 hess = zeros(noControlVar);
 for i = 1:noControlVar
  hess(i,:) = 2*J_p(:,i)'*bckGrndCovMatInv*J_p + 2*(y_p-y_o)'*bckGrndCovMatInv*H_p(:,(i-1)*noControlVar+(1:noControlVar));
 end
 hess = hess + HessPrior;

 for i = 1:noControlVar
  hess(:,i) = hess(:,i)*(ub(i)-lb(i)); hess(i,:) = hess(i,:)*(ub(i)-lb(i));
 end

end
