clear all
close all
global GelState GelSimParams
% addpath('../')


GelSimParams.Ncell = 100;
GelSimParams.dt = 0.1;
timevec = 0:GelSimParams.dt:2;
M = length(timevec)-1
GelSimParams.hx = 1/GelSimParams.Ncell;


GelState.Xcell = linspace(GelSimParams.hx/2,1-GelSimParams.hx/2,GelSimParams.Ncell)';
GelState.ThetaS = ones(GelSimParams.Ncell+2,1);

rightval = 1;
slope = 2;

C = slope*(GelState.Xcell - 1) + rightval;
plot(GelState.Xcell,C)
pause;

D = 1;
k = 2;
L = BackEulOperatorConstruct(D,GelSimParams.dt,k);
%%When it is ONLY diffusive flux that I am considering, this equates to the
%%left boundary condition Flux = k*C_x=0


for i = 1:M
time = GelSimParams.dt*i;
rhs = C;
rhs(end) = rhs(end) + rightval*D*GelSimParams.dt*2/(GelSimParams.hx^2);
Cnew = L\rhs;
plot(GelState.Xcell,Cnew)
xlim([0 1])
pause(0.1)
C = Cnew;
end
time

% true = 1 + exp(-D*pi*pi*time/4)*cos(pi*GelState.Xcell/2);
% 
% error = abs(true - C);
% 
% Linf = max(error)
% L1 = sum(error)*GelSimParams.hx
% L2 = sqrt(sum(error.^2)*GelSimParams.hx)


% rmpath('../')