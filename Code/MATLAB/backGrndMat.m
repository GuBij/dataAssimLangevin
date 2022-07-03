 function Sigma_bInv = backGrndMat(wrkDir)

 missing = -9.0;
 threshold = 0.1; % mg/m^3 Barad vol.1 p.77 (91)
 arcs=[50,100,200,400,800];
 expList = [17 18 21 22 23 24 28 29 32 37 38 39 40 41 42 46 53 54 55 56 58 59 60 70];
 noArcs = length(arcs);
 noTowers = 6;
 corrFactorsFile=readFile('PGrassEvapLossTABLE_5_4.txt',1,noArcs+1);

 noStations = 181+91*4;
 noVertStations = 9;
 noStationsTot = noStations+noTowers*noVertStations;
 innovMat = zeros(length(expList),noStationsTot);
 for expNr = expList 
  Cm = zeros(noStations,1);
  j = find( corrFactorsFile(:,1) == expNr );
  corrFactors = corrFactorsFile(j,2:end); clear j;
  filename = strcat('PPG',num2str(expNr));
  Cp = readFile(strcat(wrkDir,'conc_exp',num2str(expNr),'.txt'),0,3); Cp = Cp(:,3); 

    data=readPGFile('PGARCS.txt',0,4,true);
    j = find(data(:,1) == expNr);
    expData = data(j,:); clear data;
    for k = 1:noArcs
      j=find(expData(:,2) == k ); 
      angle = expData(j,3);
      C = corrFactors(k)*expData(j,end);
      j = find( expData(j,end) == missing ); 
      if ( length(j) > 0 )
         C(j) = NaN; clear j;
      end
      i = find(angle >= 270 & angle <= 360); angle(i) = angle(i)-360; clear i;
       if ( k == noArcs )
	innovMat(expNr,(k-1)*91 + angle + 91) = -C';
        C = Cp((k-1)*91+(1:181));
	innovMat(expNr,(k-1)*91+(1:181)) = innovMat(expNr,(k-1)*91+(1:181)) + C';
       else
	innovMat(expNr,(k-1)*91 + (angle+90)/2 + 1) = -C';
        C = Cp(((k-1)*91+1):(k*91));
	innovMat(expNr,((k-1)*91+1):(k*91)) = innovMat(expNr,((k-1)*91+1):(k*91)) + C';
       end
    end

    data=readPGFile('PGVertConc.txt',0,7,false);
    j = find(data(:,1) == expNr);
    if ( length(j) > 0 )
      expData = data(j,:); clear data j;
      for k = 1:noTowers
	C = corrFactors(2)*expData(:,k+2);
        j = find( expData(:,k+2) == missing ); 
	if ( length(j) > 0 )
	  C(j) = NaN; clear j;
	end
	C = fliplr(C');
        innovMat(expNr,noStations + (k-1)*noVertStations + (1:noVertStations)) = -C;
	C = Cp(end-noTowers*noVertStations + (k-1)*noVertStations + (1:noVertStations));
	innovMat(expNr,noStations + (k-1)*noVertStations + (1:noVertStations)) = innovMat(expNr,noStations + (k-1)*noVertStations + (1:noVertStations)) + C';
      end
    else
      innovMat(expNr,noStations + (1:(noTowers*noVertStations))) = NaN;
    end
 end

 Sigma_b = zeros(noStationsTot);
 for i=1:(noArcs-1)
  subInnovMat = innovMat(:,(i-1)*91+(1:91)); subInnovMat = subInnovMat(:);
  j = find(isnan(subInnovMat)); subInnovMat(j) = [];
  Sigma_b((i-1)*91+1,(i-1)*91+1) = mean(subInnovMat.^2);
  Sigma_b((i-1)*91*noStationsTot+(i-1)*91+(0:90)*noStationsTot+(1:91)) = Sigma_b((i-1)*91+1,(i-1)*91+1);
 end
 subInnovMat = innovMat(:,4*91+(1:181)); subInnovMat = subInnovMat(:);
 j = find(isnan(subInnovMat)); subInnovMat(j) = [];
 Sigma_b(4*91+1,4*91+1) = mean(subInnovMat.^2);
 Sigma_b(4*91*noStationsTot+4*91+(0:180)*noStationsTot+(1:181)) = Sigma_b(4*91+1,4*91+1);

 for i=1:noVertStations
  subInnovMat = innovMat(:,noStations+i+(0:(noTowers-1))*noVertStations); subInnovMat = subInnovMat(:); 
  j = find(isnan(subInnovMat)); subInnovMat(j) = [];
  Sigma_b(noStations+i,noStations+i) =  mean(subInnovMat.^2);
  Sigma_b(noStations*noStationsTot + noStationsTot*(i-1 + noVertStations*(0:(noTowers-1))) + noStations + (0:(noTowers-1))*noVertStations + i) = Sigma_b(noStations+i,noStations+i);
 end
 Sigma_bInv = Sigma_b; clear Sigma_b;
 Sigma_bInv((0:(noStationsTot-1))*noStationsTot+(1:noStationsTot)) = 1./Sigma_bInv((0:(noStationsTot-1))*noStationsTot+(1:noStationsTot));

 end
