%This function evolves a fully (electrically coupled) 
%reaction-diffusion-advection system one step using my new constrained
%Backward Euler time integrator. 
%
% function syntax:
%
%     ConstrainedBackEulStep
%
%
%     inputs:
%         none 
%     output:
%         none 


function ConstrainedBackEulStep

% Lets 'import' the two big global structs
global GelState GelSimParams

%We'll need the time step
dt = GelSimParams.dt;

%The current solvent velocity field
SolVeloc = GelState.USol;

%First, we will get ready to update the Hydrogen concentration

%The Diffusion Coefficient
D(1) = GelSimParams.Dh;
%And finally the current concentration
conccur = GelState.Hconc;

%Explicitly evaluate the advection term
advcur = AdvectionEvaluate(conccur,SolVeloc);

%Add the reaction term
explcur = GelState.HRHScur - advcur;

%Populate the entries which correspond to cells within the computational
%domain
RHSH = conccur + dt*explcur;
    
%Now we need to adjust the RHS vector to account for additional fluxes
%potentially
%%THIS CHUNK OF CODE NEEDS TO BE TESTED
RHSH(1) = RHSH(1) + dt*GelSimParams.HydFluxL/GelSimParams.hx;
RHSH(end) = RHSH(end) + GelSimParams.HydValR*D(1)*dt*2/(GelSimParams.hx^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now we'll get ready to update the Bicarbonate Concentration

%The Diffusion Coefficient
D(2) = GelSimParams.Db;
%And finally the current concentration
conccur = GelState.Bconc;

%Explicitly evaluate the advection term
advcur = AdvectionEvaluate(conccur,SolVeloc);

%Add the reaction term
explcur = GelState.BRHScur - advcur;

%Populate the entries which correspond to cells within the computational
%domain
RHSB = conccur + dt*explcur;
    
%Now we need to adjust the RHS vector to account for additional fluxes
%potentially
%%THIS CHUNK OF CODE NEEDS TO BE TESTED
RHSB(1) = RHSB(1) + dt*GelSimParams.BicFluxL/GelSimParams.hx;
RHSB(end) = RHSB(end) + GelSimParams.BicValR*D(2)*dt*2/(GelSimParams.hx^2);;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now we'll update the (negatively charged) ionic Concentration

%The Diffusion Coefficient
D(3) = GelSimParams.Di;
%And finally the current concentration
conccur = GelState.Iconc;

%Explicitly evaluate the advection term
advcur = AdvectionEvaluate(conccur,SolVeloc);

%Add the reaction term
explcur = GelState.IRHScur - advcur;

%Populate the entries which correspond to cells within the computational
%domain
RHSI = conccur + dt*explcur;

%Now we need to adjust the RHS vector to account for additional fluxes
%potentially
%%THIS CHUNK OF CODE NEEDS TO BE TESTED
RHSI(1) = RHSI(1) + dt*GelSimParams.IonFluxL/GelSimParams.hx;
RHSI(end) = RHSI(end) + GelSimParams.IonValR*D(3)*dt*2/(GelSimParams.hx^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now we'll update the (positively charged) anion Concentration

%The Diffusion Coefficient
D(4) = GelSimParams.Da;
%And finally the current and 'old' concentrations
conccur = GelState.Aconc;

%Explicitly evaluate the advection terms
advcur = AdvectionEvaluate(conccur,SolVeloc);

%Add the reaction terms
explcur = GelState.ARHScur - advcur;

%Populate the entries which correspond to cells within the computational
%domain
RHSA = conccur + dt*explcur;

%Now we need to populate the entries which correspond to ghost cells with
%the appropriate fluxes for Boundary Conditions.
RHSA(1) = RHSA(1) + dt*GelSimParams.AniFluxL/GelSimParams.hx;
RHSA(end) = RHSA(end) + GelSimParams.AniValR*D(4)*dt*2/(GelSimParams.hx^2);


val = [1,-1,1,-1];

L = ConstrainedBackEulOperatorConstruct(D,dt,val);

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