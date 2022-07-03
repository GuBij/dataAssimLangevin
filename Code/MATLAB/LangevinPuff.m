function LangevinPuff(wDir,Q);

wDir = wDir + 180;
Q = Q*1000;
%wDir = 30:0.5:60;
%Q = 56.5;

z0 = 0.008;
kappa = 0.387;
lat = 42.49333*pi/180; %58*pi/180;
radius = 100;
x0 = [0;0;0.46];
%x = [radius*cos(pi/4);radius*sin(pi/4);0.5];
uStar = 0.23; %0.39; 
L = 48; %-87;
hABL = 135; %836;

stationsFile = 'stationCoord.txt';
concFile = 'CONC_meteo5606.txt';
jacFile = 'JAC_meteo5606.txt';
hessFile = 'HESS_meteo5606.txt';

Psi = @(x) -4.7*x;
%Psi = @(x) log(0.5*(1+x.*x)*0.25.*(1+2*x+x.*x))-2*atan(x)+0.5*pi;
Psi0 = Psi(z0/L);
%Psi0 = Psi((1-15.0*z0/L)^0.25);
Umag = uStar/kappa*(log(max(x0(3),z0)/z0)-Psi(max(x0(3),z0)/L)+Psi0);
%Umag = uStar/kappa*(log(max(x0(3),z0)/z0)-Psi((1-15.0*max(x0(3),z0)/L).^0.25)+Psi0);
sigma_u = [sqrt(8.5-1.7^2)*uStar,1.7*uStar,sqrt(2.5)*uStar];
%sigma_u = [0,0.6*uStar*(-hABL/(kappa*L))^(1/3),1.4*uStar*(-x0(3)/(kappa*L))^(1/3)]; sigma_u(1) = sigma_u(2);
a = [0,0.085*sqrt(hABL*max(x0(3),z0))/sigma_u(2),0.1*hABL^0.2*(max(x0(3),z0))^0.8/sigma_u(3)]; a(1)=a(2);
%a = [0.15*hABL/sigma_u(1),0.15*hABL/sigma_u(1),0.1*x0(3)/(sigma_u(3)*(0.55+0.38*(x0(3)-z0)/L))];
a = -1./a;

data=readFile(stationsFile,0,3);
noStations = size(data,1);
C = zeros(noStations,2); C(:,1) = 0:(noStations-1);
Cjac = zeros(noStations,3); Cjac(:,1) = 0:(noStations-1);
Chess = zeros(noStations,5); Chess(:,1) = 0:(noStations-1);

T0 = 0;
T = radius/Umag;
meanDrift = @(phi) x0+Umag*(T-T0)*[cos(phi);sin(phi);0];
DispCoeff = -2./a.*sigma_u.^2.*(((exp(a*(T-T0))-2).^2-1)./a*0.5+T-T0);

SigmaInv = zeros(3);
SigmaInvJac = zeros(2);
SigmaInvHess = zeros(2);
for i = 1:noStations
   x = data(i,:); x = x';
   mu = meanDrift(wDir*pi/180);
   sinVal = sin(wDir*pi/180);
   cosVal = cos(wDir*pi/180);
   SigmaInv(1:2,1:2) = [[cosVal^2/DispCoeff(1)+sinVal^2/DispCoeff(2),sinVal*cosVal*(1/DispCoeff(1)-1/DispCoeff(2))];[sinVal*cosVal*(1/DispCoeff(1)-1/DispCoeff(2)),sinVal^2/DispCoeff(1)+cosVal^2/DispCoeff(2)]];
   SigmaInv(3,3) = 1/DispCoeff(3);
   C1 = -0.5*((x(1:2)-mu(1:2))'*SigmaInv(1:2,1:2)*(x(1:2)-mu(1:2))+(x(3)-mu(3)).^2*SigmaInv(3,3)); %sum(-0.5*(xvec-meanDrift).^2./var);
   C2 = -0.5*((x(1:2)-mu(1:2))'*SigmaInv(1:2,1:2)*(x(1:2)-mu(1:2))+(x(3)+mu(3)).^2*SigmaInv(3,3));
   C(i,2) = Q*(exp(C1)+exp(C2))/sqrt(8*pi^3*prod(DispCoeff)); %prod(sqrt(2*pi*var));

   muJac = Umag*(T-T0)*[-sinVal;cosVal];
   muHess = Umag*(T-T0)*[-cosVal;-sinVal];
   SigmaInvJac = [[2*sinVal*cosVal*(-1/DispCoeff(1)+1/DispCoeff(2)),(2*cosVal^2-1)*(1/DispCoeff(1)-1/DispCoeff(2))];[(2*cosVal^2-1)*(1/DispCoeff(1)-1/DispCoeff(2)),2*sinVal*cosVal*(1/DispCoeff(1)-1/DispCoeff(2))]];
   SigmaInvHess = [[2*(2*cosVal^2-1)*(-1/DispCoeff(1)+1/DispCoeff(2)),-2*(2*sinVal*cosVal)*(1/DispCoeff(1)-1/DispCoeff(2))];[-2*(2*sinVal*cosVal)*(1/DispCoeff(1)-1/DispCoeff(2)),2*(2*cosVal^2-1)*(1/DispCoeff(1)-1/DispCoeff(2))]];
   Cjac(i,2) = C(i,2)*(muJac'*SigmaInv(1:2,1:2)*(x(1:2)-mu(1:2))-0.5*(x(1:2)-mu(1:2))'*SigmaInvJac*(x(1:2)-mu(1:2)));
   Chess(i,2) = C(i,2)*((muJac'*SigmaInv(1:2,1:2)*(x(1:2)-mu(1:2))-0.5*(x(1:2)-mu(1:2))'*SigmaInvJac*(x(1:2)-mu(1:2)))^2-muJac'*SigmaInv(1:2,1:2)*muJac+muHess'*SigmaInv(1:2,1:2)*(x(1:2)-mu(1:2))+2*muJac'*SigmaInvJac*(x(1:2)-mu(1:2))-0.5*(x(1:2)-mu(1:2))'*SigmaInvHess*(x(1:2)-mu(1:2)));
end
Cjac(:,2:end) = Cjac(:,2:end)*pi/180;
Chess(:,2:end) = Chess(:,2:end)*(pi/180)^2;

fileID=fopen(sprintf('output/%s',concFile),'w');
fprintf(fileID,'\n');
fprintf(fileID,'0\t\t%d\t\t%6.5e\n',C');
fclose(fileID);

fileID=fopen(sprintf('output/%s',jacFile),'w');
fprintf(fileID,'\n');
fprintf(fileID,'0\t\t%d\t\t%6.5e\t\t%6.5e\n',Cjac');
fclose(fileID);

fileID=fopen(sprintf('output/%s',hessFile),'w');
fprintf(fileID,'\n');
fprintf(fileID,'0\t\t%d\t\t%6.5e\t\t%6.5e\t\t%6.5e\t\t%6.5e\n',Chess');
fclose(fileID);

%{
figure(1);
plot(wDir,C,'-k'); hold off;
ylabel('$c$','FontSize',18,'interpreter','latex');
figure(2);
plot(wDir,Cjac,'-k'); hold off;
ylabel('$\partial c/\partial \phi$','FontSize',18,'interpreter','latex');
figure(3);
plot(wDir,Chess,'-k'); hold off;
ylabel('$\partial ^{2} c/\partial \phi ^{2}$','FontSize',18,'interpreter','latex');

xlabel('$\phi$ [\circ]','interpreter','latex','FontSize',18);
ylabel('$c$ [mg/m$^3$]','interpreter','latex','FontSize',18);
title(sprintf('%d $m$ - arc',arcs(k)),'FontSize',18,'interpreter','latex');
set(gcf,'renderer','painters','Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(gcf,strcat(filename,'_arc',num2str(k)),'-dpdf'); %'-dpng');
%}

end

