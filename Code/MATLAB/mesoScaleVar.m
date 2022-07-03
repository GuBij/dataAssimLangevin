function [sigma_meso,turbVarStdDev] = mesoScaleVar(sigmaAzim,modelParams,propConstStable)

  LCol = 2;
  uStarCol = 3;
  hCol = 4;
  sigma_uCol = 6;
  sigma_vCol = 7;
  UCol = 10;
  lat = 42.49333; 
  kappa = 0.387;
  kappaCorr = round(kappa/0.35,2);
  missing = -999.0;

  if ( nargin == 2 )
    propConstStable = 1.7;
  end

  modelParams(uStarCol) = kappaCorr*modelParams(uStarCol);
  if ( modelParams(hCol) == missing && abs(modelParams(LCol))>=500 )
    fCori = 2*7.2921e-5*sin(lat*pi/180);
    modelParams(hCol) = 0.1*modelParams(uStarCol)/fCori;
  elseif ( modelParams(hCol) == missing && modelParams(LCol) > 0 && modelParams(LCol) < 500 )
    fCori = 2*7.2921e-5*sin(lat*pi/180);
    modelParams(hCol) = 0.4*sqrt(modelParams(LCol)*modelParams(uStarCol)/fCori);
  elseif (  modelParams(hCol) == missing )
    error('Pasquill did not estimate mixingheight for experiment %d.',expNr);
  end

  if ( abs(modelParams(LCol))>= 500 )
    modelParams(sigma_vCol) = sqrt(4.1)*modelParams(uStarCol);
  elseif ( modelParams(LCol) < 0 )
    modelParams(sigma_vCol) = 0.6*modelParams(uStarCol)*(modelParams(hCol)/(kappa*(-modelParams(LCol))))^(1/3);;
  elseif ( modelParams(LCol) > 0 )
    modelParams(sigma_uCol) = sqrt(8.5-propConstStable^2)*modelParams(uStarCol);
    modelParams(sigma_vCol) = propConstStable*modelParams(uStarCol);
  end
  turbVarStdDev = modelParams(sigma_vCol)/modelParams(UCol)*180/pi;
  sigma_meso = sigmaAzim^2-turbVarStdDev^2;
  if ( sigma_meso > 0 )
	sigma_meso = sqrt(sigma_meso);
  else
	sigma_meso = 0;
	warning('Mesoscale variability is not present anymore!');
  end

end
