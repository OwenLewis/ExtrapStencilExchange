clear all
close all
global GelState GelSimParams

refine = 2;

GelSimParams.Ncell = 100*refine;
GelSimParams.Nedges = GelSimParams.Ncell-1;
GelSimParams.dt = 0.00001;
timevec = 0:GelSimParams.dt:1;
M = length(timevec)-1
GelSimParams.hx = 1/GelSimParams.Ncell;


GelState.Xcell = linspace(GelSimParams.hx/2,1-GelSimParams.hx/2,GelSimParams.Ncell)';
GelState.XcellExtend = linspace(-GelSimParams.hx/2,1+GelSimParams.hx/2,GelSimParams.Ncell+2)';
GelState.Xedge = linspace(GelSimParams.hx,1-GelSimParams.hx,GelSimParams.Nedges)';
GelState.XedgeExtend = linspace(0,1,GelSimParams.Nedges+2)';
C = 1 + cos(pi.*GelState.Xcell/2);
% C = ones(size(GelState.Xcell));



GelState.ThetaS = ones(GelSimParams.Ncell+2,1);
GelState.Hconc = C;
GelState.Hold = C;
GelState.Bconc = C;
GelState.Bold = C;
GelState.Iconc = 0*C;
GelState.Iold = 0*C;
GelState.Aconc = 0*C;
GelState.Aold = 0*C;
GelState.DPsi = 0*GelState.Xedge;

GelSimParams.Kbind = 0;

GelSimParams.Dh = 0.1;
GelSimParams.Db = 0.2;
GelSimParams.Da = 0;
GelSimParams.Di = 0;


GelSimParams.HydExchangeRate = 0;
GelSimParams.HydExchangerParam = 0;
GelSimParams.BicExchangeRate = 0;
GelSimParams.BicExchangerParam = 0;

HydRight = 1;
BicRight = 1;

GelSimParams.HydValR = HydRight;
GelSimParams.BicValR = BicRight;
GelSimParams.IonValR = 0;
GelSimParams.AniValR = 0;


% plot(GelState.Xcell,C)
% pause;


% rhs(end) + rightval*D*GelSimParams.dt*2/(GelSimParams.hx^2);

disp('starting time loop')
tic
for i = 1:M
time = GelSimParams.dt*i;

RHSH = GelState.Hconc;
RHSH(end) = RHSH(end) + HydRight*GelSimParams.Dh*GelSimParams.dt*2/(GelSimParams.hx^2);

RHSB = GelState.Bconc;
RHSB(end) = RHSB(end) + BicRight*GelSimParams.Db*GelSimParams.dt*2/(GelSimParams.hx^2);

RHSI = GelState.Iconc;
RHSI(end) = RHSI(end) + GelSimParams.IonValR*GelSimParams.Di*GelSimParams.dt*2/(GelSimParams.hx^2);

RHSA = GelState.Aconc;
RHSA(end) = RHSA(end) + GelSimParams.AniValR*GelSimParams.Da*GelSimParams.dt*2/(GelSimParams.hx^2);

val = [1,-1,1,-1];
D = [GelSimParams.Dh,GelSimParams.Db,GelSimParams.Di,GelSimParams.Da];

L = ConstrainedBackEulOperatorConstruct(D,GelSimParams.dt,val);

RHS = [RHSH;RHSB;RHSI;RHSA;zeros(GelSimParams.Nedges,1)];


newconcs = L\RHS;

%Lets pluck out the entries corresponding to hydrogen
concnew = newconcs(1:GelSimParams.Ncell);
%Current becomes old, and new becomes current
GelState.Hold = GelState.Hconc;
GelState.Hconc = concnew;

%Now we'll pluck out the entries corresponding to bicarbonate
concnew = newconcs(GelSimParams.Ncell+1:2*GelSimParams.Ncell);
%Current becomes old, and new becomes current
GelState.Bold = GelState.Bconc;
GelState.Bconc = concnew;

%Now we'll pluck out the entries corresponding to negative Ions
concnew = newconcs(2*GelSimParams.Ncell+1:3*GelSimParams.Ncell);
%Current becomes old, and new becomes current
GelState.Iold = GelState.Iconc;
GelState.Iconc = concnew;

%Now we'll pluck out the entries corresponding to positive Anions
concnew = newconcs(3*GelSimParams.Ncell+1:4*GelSimParams.Ncell);
%Current becomes old, and new becomes current
GelState.Aold = GelState.Aconc;
GelState.Aconc = concnew;

%And finally, lets pluck out the electric potential gradient
GelState.DPsi = newconcs(4*GelSimParams.Ncell+1:end);
end
toc
time;

Deff = 2*GelSimParams.Dh*GelSimParams.Db/(GelSimParams.Dh + GelSimParams.Db);

true = 1 + exp(-Deff*pi*pi*time/4)*cos(pi*GelState.Xcell/2);


truePhi = ((GelSimParams.Db - GelSimParams.Dh)/(GelSimParams.Db + GelSimParams.Dh))*(-pi*exp(-pi*pi*Deff*time/4).*sin(pi*GelState.Xedge/2)/2)./(1 +exp(-pi*pi*Deff*time/4).*cos(pi*GelState.Xedge/2));

errorH = true - GelState.Hconc;
disp('Hydrogen')
Linf = max(abs(errorH))
L1 = sum(abs(errorH))*GelSimParams.hx
L2 = sqrt(sum(errorH.^2)*GelSimParams.hx)

errorPhi = truePhi - GelState.DPsi;
disp('Grad Psi')
Linf = max(abs(errorPhi))
L1 = sum(abs(errorPhi))*GelSimParams.hx
L2 = sqrt(sum(errorPhi.^2)*GelSimParams.hx)




figure(1)
plot(GelState.Xcell,GelState.Hconc,'--',GelState.Xcell,GelState.Bconc,'--',GelState.Xcell,GelState.Aconc,'--',GelState.Xcell,GelState.Iconc,'--',GelState.Xcell,true,'LineWidth',3)
xlabel('X')
ylabel('Concentration')
title('Time t = 2')
legend('Hyd.','Bicarbonate','Chloride','Sodium','True Soln.','location','southwest')
set(gca,'FontSize',14)

figure(2)
plot(GelState.Xedge,GelState.DPsi,'--',GelState.Xedge,truePhi,'LineWidth',3)
xlabel('X')
ylabel('Grad Psi')
title('Time t = 2')
legend('Computed','True Soln.','location','southwest')
set(gca,'FontSize',14)
