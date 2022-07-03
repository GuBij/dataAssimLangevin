% run function like this: bayOpt(0,'windDirection',170,'sourceStrength',50);

function bayOpt(expNr,varargin) 

 global bckGrndCovMatInv priorWeights y_o x_o sigmaAzimMin varPresent concFile sensFile hessFile inputFileU inputFileK inputFile rmStations modelParams inputFileUOrig inputFileKOrig 

 observations='predictions'; % 'measurements'; 'predictions' 
 maxNoControlVar = 3;
 error_wDir = 10;
 pct = 0.05;

 if isstr(expNr)
  expNr = str2num(expNr);
 end

 if ( nargin < 3 | rem(nargin-1,2) ~= 0 )
   error('Not enough input arguments: provide experiment number and pairs of control parameter names and values.');
 end
 noControlVar = (nargin-1)/2;
 x_o = zeros(noControlVar,1);
 lb = zeros(noControlVar,1);
 ub = zeros(noControlVar,1);
 priorWeights = zeros(1,noControlVar);
 varID = zeros(noControlVar,1);
 varPresent = zeros(1,maxNoControlVar);
 modelParams = 0;
 sigmaAzimMin = 0;
 for i = 1:noControlVar
   if isstr(varargin{2*i})
     varargin{2*i} = str2num(varargin{2*i});
   end

   if strcmp(varargin{2*i-1},'windDirection')
     x_o(i) = varargin{2*i};
     lb(i) = -Inf; 
     ub(i) = Inf; 
     varID(i) = 1;
     varPresent(1) = 1;
   elseif strcmp(varargin{2*i-1},'windDirStdDev')
     x_o(i) = varargin{2*i};
     modelParams=readFile('modelParamList.txt',0,10);
     expFileId = find( modelParams(:,1) == expNr );
     if ( length(expFileId) == 0 )
       error('experiment %d does not occur in the parameter list.',expNr);
     else
       modelParams = modelParams(expFileId,:);
     end
     clear expFileId;
     [sigmaAzimMeso,sigmaAzimMin] = mesoScaleVar(x_o(i),modelParams);
     stdDevAzim = x_o(i); 
     if ( sigmaAzimMeso == 0 )
       ub(i) = 12.5; %7.5; %Inf;
       x_o(i) = 0.5*(sigmaAzimMin+ub(i));
       disp(sprintf('\nInitial value for windDirStdDev has been changed to %4.2f\n',x_o(i)));
     else
       ub(i) = Inf;
       priorWeights(i) = 2/(x_o(i)-sigmaAzimMin);
     end
     lb(i) = sigmaAzimMin;
     varID(i) = 2;
     varPresent(2) = 1;
   elseif strcmp(varargin{2*i-1},'sourceStrength')
     x_o(i) = varargin{2*i};
     ub(i) = Inf;
     varID(i) = 3;
     priorWeights(i) = 2/x_o(i);
     varPresent(3) = 1;
   else
      error('Variable "%s" is not included. Included variables: "windDirection", "windDirStdDev" and "sourceStrength".',varargin{2*i-1});
   end
 end
 [varID,orderVar]=sort(varID); 
 x_o = x_o(orderVar);
 lb = lb(orderVar);
 ub = ub(orderVar);
 priorWeights = priorWeights(orderVar); clear orderVar;
 varPresent = logical(varPresent);
 if ( varPresent(1) & varPresent(2) )
   priorWeights(1) = 1/stdDevAzim^2;
 elseif varPresent(1)
   priorWeights(1) = (norminv(1-pct)/error_wDir)^2;
 end

 concDir = 'output/PPG_stable/';
 concFile = 'output/CONC_meteo5606.txt';
 sensFile = 'output/JAC_meteo5606.txt';
 hessFile = 'output/HESS_meteo5606.txt';
 inputFileUOrig = sprintf('PARAMS_UopenField_meteo5606_exp%s.txt',num2str(expNr)); 
 inputFileU = 'PARAMS_UopenField_meteo5606.txt'; 
 inputFileKOrig = sprintf('PARAMS_KTaylorOF_meteo5606_exp%s.txt',num2str(expNr));
 inputFileK = 'PARAMS_KTaylorOF_meteo5606.txt';
 inputFile = 'inputDict';
 obsFile = 'output/conc_'; 
 bckGrndCovMatInv = backGrndMat(concDir);

 if strcmp(observations,'predictions')
   y_o = readFile(strcat(obsFile,'exp',num2str(expNr),'.txt'),0,3); y_o = y_o(:,3); 
   rmStations=[];
 elseif strcmp(observations,'measurements')
   noArcs = 5;
   noTowers = 6;
   noVertStations = 9;
   corrFactors=readFile('PGrassEvapLossTABLE_5_4.txt',1,noArcs+1);
   j = find( corrFactors(:,1) == expNr );
   corrFactors = corrFactors(j,2:end);
   arcData=readPGFile('PGARCS.txt',0,4,true);
   j = find(arcData(:,1) == expNr);
   expData = arcData(j,:); clear arcData;
   y_o = zeros(181+91*4+noTowers*noVertStations,1);
   rmStations = [];
   for k = 1:noArcs
     j=find(expData(:,2) == k );
     angle = expData(j,3);
     i = find(angle >= 270 & angle <= 360); angle(i) = angle(i)-360;
     i = find( expData(j,end) == missing ); 
     if ( k == noArcs )
       jj = (k-1)*91 + angle + 91;
       ii = (k-1)*91 + angle(i) + 91;
     else
       jj = (k-1)*91 + (angle+90)/2 + 1;
       ii = (k-1)*91 + (angle(i)+90)/2 + 1;
     end
     y_o(jj) = corrFactors(k)*expData(j,end); clear j jj;
     if ( length(ii) > 0 )
      rmStations = [rmStations;ii]; clear i ii;
     end
   end

   towerData=readPGFile('PGVertConc.txt',0,noTowers+1,false);
   j = find(towerData(:,1) == expNr);
   if ( noVertStations > 0 & length(j) > 0 )
     expData = towerData(j,:); clear towerData j;
     for k = 1:noTowers
       C = corrFactors(2)*expData(:,k+2);
       j = find( expData(:,k+2) == missing );
       C = fliplr(C');
       if ( length(j) > 0 )
         rmStations = [rmStations;181+91*4+k*noVertStations+1-j];
       end
       y_o(181+91*4+(k-1)*noVertStations+(1:noVertStations)') = C';
     end
   else
     rmStations = [rmStations;181+91*4+(1:(noVertStations*noTowers))'];
   end
   clear expData;
   y_o(rmStations) = []; bckGrndCovMatInv(rmStations,:) = []; bckGrndCovMatInv(:,rmStations) = [];
 end

 options = optimoptions('fmincon','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'HessianFcn','objective','MaxIterations',200,'Display','iter','MaxFunctionEvaluations',10^4); 
 [x,fval,exitflag,output] = fmincon(@objFuncBay,x_o,[],[],[],[],lb,ub,[],options); 

 disp(sprintf('\nExperiment %d:\n',expNr));
 disp(output);
 disp(exitflag);

 for i = 1:noControlVar
   if ( varID(i) == 1 )
     disp(sprintf('\nThe wind direction is %5.2f degrees.',x(i)));
   elseif ( varID(i) == 2 & stdDevAzim < sigmaAzimMin & round(x(i),2) == round(sigmaAzimMin,2) )
     disp(sprintf('\nThe standard deviation of the wind direction is %4.2f degrees (%4.2f degrees).',stdDevAzim,sigmaAzimMin));
   elseif ( varID(i) == 2 )
     disp(sprintf('\nThe standard deviation of the wind direction is %5.3f degrees.',x(i)));
   elseif ( varID(i) == 3 )
     disp(sprintf('\nThe source strength is %5.2f g/s.',x(i)));
   end
 end
 disp(sprintf('\n'));

end 
