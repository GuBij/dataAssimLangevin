function [FB,FBfp,FBfn,NMSE,FAC2,FAC10,n]=predSkill()

 threshold = 0.1; % mg/m^3 Barad vol.1 p.77 (91)
 caseList = [17,18,21,22,23,24,28,29,32,37,38,39,40,41,42,46,53,54,55,56,58,59,60,70];
% caseList = [17,18,21,22,23,24,28,29,32,37,38,39,40,41,42,46,53,54,55,56,58,59,60,70]; %stable cases
% caseList = [5,6,7,8,9,10,11,12,15,16,19,20,25,26,27,30,31,33,34,43,44,45,48,49,50,51,57,61,62]; %unstable cases
 wrkDir = 'PPG_stable/'; 
 Cm = [];
 Cp = [];
 noCases = length(caseList);
 MR = zeros(noCases*5,1);
 for k = 1:noCases
  Cp_ = readFile(strcat(wrkDir,'conc_exp',num2str(caseList(k)),'.txt'),0,3); Cp_ = Cp_(:,3); %Cp_ = Cp_(1:(end-54),3);
  Cp = [Cp;Cp_];

  Cm_ = PGConc(caseList(k),true,{wrkDir}); 
  Cm = [Cm;Cm_];
  Cm_ = PGConc(caseList(k),false,{wrkDir});
  if ( length(Cm_) == 0 )
    Cp((end-53):end) = [];
  end
  Cm = [Cm;Cm_]; 
 end
 clear Cm_ Cp_; 
 tooLow=find( Cp < threshold & Cm < threshold );
 Cp(tooLow,:) = [];
 Cm(tooLow,:) = [];

 n = length(Cm);
 sumCp = sum(Cp);
 sumCm = sum(Cm);
 ifp = find(Cp>Cm);
 ifn = find(Cm>Cp);
 FB = 2*(sumCm - sumCp)/(sumCm+sumCp);
 FBfp = (sum(abs(Cp(ifp)-Cm(ifp)) + Cp(ifp)-Cm(ifp)))/(sumCm+sumCp);
 FBfn = (sum(abs(Cp(ifn)-Cm(ifn)) + Cm(ifn)-Cp(ifn)))/(sumCm+sumCp);
 NMSE = n*sum((Cm-Cp).^2)/(sumCm*sumCp);
 ratio = Cp./Cm;
 FAC2 = length(find(ratio >= 0.5 & ratio <= 2.0 ))/n;
 FAC10 = length(find(ratio >= 0.1 & ratio <= 10.0 ))/n;
 NMSEs = 4*FB^2/(4-FB^2);

end
