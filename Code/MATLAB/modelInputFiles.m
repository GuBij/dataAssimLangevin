function modelInputFiles(expNr,heff)

  if isstr(expNr)
    expNr = str2num(expNr);
  end

  if isstr(heff)
    heff = str2num(heff);
  end

  expNrCol = 1;
  LCol = 2;
  uStarCol = 3;
  hCol = 4;
  wDirCol = 5;
  sigma_uCol = 6;
  sigma_vCol = 7;
  sigma_wCol = 8;
  sigma_wDirCol = 9;
  UCol = 10;
  zref = 2;
  lat = 42.49333; 
  kappa = 0.387;
  kappaCorr = round(kappa/0.35,2);
  missing = -999.0;

  modelParams=readFile('modelParamList.txt',0,10);
  expFileId = find( modelParams(:,1) == expNr );
  if ( length(expFileId) == 0 )
    error('experiment %d does not occur in the parameter list.',expNr);
  else
    modelParams = modelParams(expFileId,:);  
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

  propConst = 0;
  if ( abs(modelParams(LCol))>= 500 )
    modelParams(sigma_uCol) = sqrt(6.3)*modelParams(uStarCol);
    modelParams(sigma_vCol) = sqrt(4.1)*modelParams(uStarCol);
    modelParams(sigma_wCol) = sqrt(1.7)*modelParams(uStarCol);
  elseif ( modelParams(LCol) < 0 )
    propConst = 1.4*modelParams(uStarCol)/(kappa*(-modelParams(LCol)))^(1/3); %modelParams(sigma_wCol)/zref^(1/3);
    modelParams(sigma_uCol) = 0.6*modelParams(uStarCol)*(modelParams(hCol)/(kappa*(-modelParams(LCol))))^(1/3);
    modelParams(sigma_vCol) = modelParams(sigma_uCol);
  elseif ( modelParams(LCol) > 0 )
    modelParams(sigma_uCol) = sqrt(5.6)*modelParams(uStarCol);
    modelParams(sigma_vCol) = 1.7*modelParams(uStarCol);
    modelParams(sigma_wCol) = sqrt(2.5)*modelParams(uStarCol);
  end
  sigma_meso = modelParams(sigma_wDirCol)^2-(modelParams(sigma_vCol)/modelParams(UCol)*180/pi)^2;
  if ( sigma_meso > 0 )
	sigma_meso = sqrt(sigma_meso);
	warning('mesoscale variability of %3.1f degrees need to be added for exp. %d',sigma_meso,expNr);
  elseif ( -sigma_meso > 1 )
	warning('estimated turbulent variability deviates with %3.1f degrees for exp. %d',sqrt(-sigma_meso),expNr);
	sigma_meso = 0;
  else
	sigma_meso = 0;
  end

  fid1 = fopen('output/PARAMS_UopenField_meteo5606.txt','w');
  fprintf(fid1,'0\t%3.2f\t%5.1f\t%5.2f\t%d\t%3.1f\n',modelParams(uStarCol),modelParams(LCol),heff,modelParams(wDirCol),sigma_meso);
  fclose(fid1);

  fid1 = fopen('output/PARAMS_KTaylorOF_meteo5606.txt','w');
  fprintf(fid1,'0\t%5.4f\t%5.4f\t%5.4f\t%4.2f\t%3.2f\t%5.1f\t%5.2f\t%d',modelParams(sigma_uCol),modelParams(sigma_vCol),modelParams(sigma_wCol),propConst,modelParams(uStarCol),modelParams(LCol),heff,round(modelParams(hCol)));
  fclose(fid1);

end
