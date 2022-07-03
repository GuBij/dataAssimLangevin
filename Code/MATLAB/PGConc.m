%run function like this: PGConc(41,true,{'PPG_stable/','Bayopt/'},[0,35]); or PGConc(41,false,{'PPG_stable/','Bayopt/'});
%arcBool = true plots the results for the arcs from experiment expNr, arcBool = false plots the results for the arcs
%wrkDir = working directory
%xRange = fixes the ranges of the x-axis over the different plots
%fineResolution = 'stable' or 'unstable' is used for plotting the results of simulations using a finer mesh

function Cm = PGConc(expNr,arcBool,wrkDir,xRange,fineResolution)

 missing = -9.0;
 threshold = 0.1; % mg/m^3 Barad vol.1 p.77 (91)
 boolPredBar = false; % error bar on predictions?
 arcs=[50,100,200,400,800];
 noArcs = length(arcs);
 corrFactors=readFile('PGrassEvapLossTABLE_5_4.txt',1,noArcs+1);
 j = find( corrFactors(:,1) == expNr );
 corrFactors = corrFactors(j,2:end); clear j;
 filename = strcat('thesisPlots/PPG',num2str(expNr)); % directory where to save the picture
 lineStyle={'--','-'};
 Cp_mat = [];
 for i = 1:length(wrkDir)
  Cp = readFile(strcat(wrkDir{i},'conc_exp',num2str(expNr),'.txt'),0,3); Cp = Cp(:,3); % 'CONC1e04_meteo560814_fine.txt',0,3); Cp = Cp(:,3);
  Cp_mat = [Cp_mat,Cp];
 end

 if boolPredBar
%  wrkDirPredBar = {'temp/Table4/stable/KS/np1e3_dt0.05/','temp/Table4/stable/KS/np5e4_dt0.025/','temp/Table4/stable/PI/np1e3_dt0.05/','temp/Table4/stable/PI/np5e4_dt0.025/'};
  wrkDirPredBar = {'temp/Table4/unstable/KS/np5e4_dt0.005/','temp/Table4/unstable/PI/np5e4_dt0.005/'};
  predBar = [];
  lineStylePredBar={'.r','.k'};
  markerStyle={'+','*'};
  for i = 1:length(wrkDirPredBar)
   Cp = readFile(strcat(wrkDirPredBar{i},'conc_exp',num2str(expNr),'.txt'),0,3); Cp = Cp(:,3);
   predBar = [predBar,Cp];
  end
  clear wrkDirPredBar;
 end

 clear Cp;
 if ( arcBool )
   Cm = zeros(181+91*4,1); 
 else
   Cm = zeros(6*9,1);
 end

 if ( arcBool )
   data=readPGFile('PGARCS.txt',0,4,true);
   j = find(data(:,1) == expNr);
   expData = data(j,:); clear data;
   for k = 1:noArcs
     j=find(expData(:,2) == k );
     angle = expData(j,3);
     C = corrFactors(k)*expData(j,end);
     j = find( expData(j,end) == missing | expData(j,end) < threshold );
     if ( length(j) > 0 )
        angle(j) = []; C(j) = []; clear j;
     end
     i = find(angle >= 270 & angle <= 360); angle(i) = angle(i)-360; clear i;
     if ( length(angle) > 0 )
      figure(k);
      errorbar(angle,C,0.1*C,'ok','MarkerSize',12,'LineWidth',1.1); hold on;
      if ( k == noArcs )
       Cm((k-1)*91 + angle + 91) = C;
       if ( nargin == 5 & strcmp(fineResolution,'stable'))
	elemPrevArc = 101;
	angle = (-30):0.6:30; % 10:0.6:70; %(-30):0.6:30;
       else
	elemPrevArc = 91;
	angle = (-90):90;
       end
       for l = 1:length(wrkDir)
        C = Cp_mat((k-1)*elemPrevArc+(1:length(angle)),l); %Cp((k-1)*101+(1:101)); %Cp((k-1)*91+(1:181));
        i = find( C < threshold ); C(i) = []; angle_cpy = angle; angle_cpy(i) = [];
        plot(angle_cpy,C,strcat(lineStyle{l},'k'),'LineWidth',1.1); hold on;
	if boolPredBar
	  C = predBar((k-1)*91+(1:181),l); %3*(l-1)+(1:3));
	  i = find( C(:,1) < threshold ); C(i,:) = []; angle_cpy = (-90):90; angle_cpy(i) = [];
	  plot(angle_cpy,C,strcat(markerStyle{l},'k'),'MarkerSize',12); hold on;
	end

       end
       hold off; 
       clear angle_cpy angle C i;
      else
       Cm((k-1)*91 + (angle+90)/2 + 1) = C;
       if ( nargin == 5 & strcmp(fineResolution,'stable'))
        angle = (-30):0.6:30; % 10:0.6:70; %(-30):0.6:30;
       else
        angle = (-90):2:90;
       end
       for l = 1:length(wrkDir)
        C = Cp_mat(((k-1)*length(angle)+1):(k*length(angle)),l); %Cp((k-1)*101+(1:101)); %Cp(((k-1)*91+1):(k*91));
        i = find( C < threshold ); C(i) = []; angle_cpy = angle; angle_cpy(i) = [];
        plot(angle_cpy,C,strcat(lineStyle{l},'k'),'LineWidth',1.1); hold on;

        if boolPredBar
	  C = predBar(((k-1)*91+1):(k*91),l); %3*(l-1)+(1:3));
	  i = find( C(:,1) < threshold ); C(i,:) = []; angle_cpy = (-90):2:90; angle_cpy(i) = [];
	  plot(angle_cpy,C,strcat(markerStyle{l},'k'),'MarkerSize',12); hold on;
        end
       end
       hold off;
       clear angle_cpy angle C i;
      end
      if ( nargin >= 4 & length(xRange) == 2 )
        xlim(xRange);
      end
      set(gca,'FontSize',20);
      xlabel('angle $[^{\circ}]$','interpreter','latex','FontSize',30);
      ylabel('$c$ [mg/m$^3$]','interpreter','latex','FontSize',30);
 
     set(gcf,'renderer','painters','Units','Inches');
     pos = get(gcf,'Position');
     set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%     print(gcf,strcat(filename,'_arc',num2str(k)),'-dpdf'); % '-dpng'); % '-dpdf');
     end
   end
 else 
   noTowers = 6;
   data=readPGFile('PGVertConc.txt',0,7,false);
   j = find(data(:,1) == expNr);
   if ( length(j) == 0 )
     warning('Tower measurements were not conducted for experiment %d.',expNr);
     Cm = []; return
   end
   expData = data(j,:); clear data j;
   z = expData(:,2); 
   if ( nargin == 5 )
     zp = 0.25:0.25:17.5;
   else
     zp = z; zp = fliplr(zp');
   end
   for k = 1:noTowers
     j = find( expData(:,k+2) == missing | expData(:,k+2) < threshold );
     C = corrFactors(2)*expData(:,k+2);
     Cm(((k-1)*9+1):(k*9)) = C;
     z_cpy = z;
     if ( length(j) > 0 )
        z_cpy(j) = []; C(j) = []; Cm((k-1)*9 + j) = 0; clear j; 
     end
     if ( length(C) > 0 )
      figure();
      for l = 1:length(wrkDir)
       plot(Cp_mat(end-noTowers*length(zp)+(((k-1)*length(zp)+1):(k*length(zp))),l),zp,strcat(lineStyle{l},'k'),'LineWidth',1.1); hold on;
       if boolPredBar
	zp2 = z; zp2 = fliplr(zp2');
	plot(predBar(end-noTowers*length(zp2)+(((k-1)*length(zp2)+1):(k*length(zp2))),l),zp2,strcat(markerStyle{l},'k'),'MarkerSize',12); hold on;
       end
      end
      errorbar(C,z_cpy,0.1*C,'horizontal','ok','MarkerSize',12,'LineWidth',1.1); hold off;
      set(gca,'FontSize',20);
      ylabel('height [m]','interpreter','latex','FontSize',30);
      xlabel('$c$ [mg/m$^3$]','interpreter','latex','FontSize',30);
%      title(sprintf('Exp. %d, tower %d',expNr,k),'FontSize',18,'interpreter','latex');
%      legend({'without DA','with DA','measurements'},'Location','northeast');
     
      set(gcf,'renderer','painters','Units','Inches');
      pos = get(gcf,'Position');
      set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%      print(gcf,strcat(filename,'_tower',num2str(k)),'-dpng'); % '-dpng'); % '-dpdf');
     end
   end
 end

end
