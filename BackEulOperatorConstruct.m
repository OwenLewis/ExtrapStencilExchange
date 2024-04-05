%This function construts the ``gell averaged'' operators for Backwards
%Euler iteration of diffusive species in the two-phase gel model. These 
%operators are then used to advance the reaction-diffusion-advection system. 
%For the purposes of this construction, we assume that the R-D-A system is 
%subject to Robin boundary condtions at the left and homogeneous Neumann 
%conditions at the right. We also assume that the solvent volume fraction 
%is subject to a given advective flux at the left and outflow boundary 
%conditions at the right.
%
% function syntax:
%
%     [Lold,Lcur,Lfut] = BackEulOperatorConstruct(DiffCoeff,dt,BndFluxCoeff)
%
%
%     inputs:
%         DiffCoeff is a number which defines the diffusion coefficient for
%           the species which these operators will update the heat equation
%           for
%         dt is the time update step size. 
%         BndFluxCoeff is a scalar which represents the coefficient which
%           appears in front of the wall concentration of this species in
%           the boundary flux calculation
%     output:
%         Lfut is the Ncell x Ncell sparce matrix which will be
%           applies  the operation of the SBD scheme to the future time
%           step. This one will be inverted to solve for the new values.


function Lfut = BackEulOperatorConstruct(DiffCoeff,dt,BndFluxCoeff)

%Lets 'import' the two big global structs
global GelState GelSimParams
% 'time to form gel velocity operator'

%Here are some parameters we need to define the sizes of things
hx = GelSimParams.hx;
Ncell = GelSimParams.Ncell;

%And the diffusion coefficient, whos name i don't want to type over & over
D = DiffCoeff;



%IMPORTANT, WE ASSUME THAT GELSTATE.THETAS IS ALREADY OF SIZE NEDGE+2, AND
%THEREFORE CONTAINS THE GHOST CELLS WHICH HAVE BEEN CALCULATED TO REFLECT
%THE BOUNDARY CONDITIONS ON SOLVENT VOLUME FRACTION

%For ease of typing, we will make a locally scoped array for volume
%fraction
workingSol = GelState.ThetaS;
edgeSol = (GelState.ThetaS(1:end-1) + GelState.ThetaS(2:end))/2;
edgeSol(1) = 0;

%Puting together the Laplacian-like operator for the Solvent
LSolUDiag = D*(edgeSol(2:end-1))./(workingSol(2:end-2)*hx^2);
LSolLDiag = D*(edgeSol(2:end-1))./(workingSol(3:end-1)*hx^2);
LSolMDiag = -D*(edgeSol(1:end-1) + edgeSol(2:end))./(workingSol(2:end-1)*hx^2);

%This is an adjustment to the last entry of the main diagonal to account
%for dirichlet BC on a cell centered grid
LSolMDiag(end) = LSolMDiag(end) -D*edgeSol(end)/(workingSol(end-1)*hx^2);

%Now we need to adjust the first row to account for flux through the left
%boundary
factor = BndFluxCoeff*mean(workingSol(1:2))/(workingSol(2)*hx);

LSolUDiag(1) = LSolUDiag(1) - factor/2;
LSolMDiag(1) = LSolMDiag(1) + 3*factor/2;
%%%The above adjustment basically specifies that there is a flux at the
%%%left boundary which is given by F = BndFluxCoeff*C|x=0. We are
%%%approximating C|x=0 via linear extrapolation from the first two interior
%%%cell centers. 


% keyboard
%Now we will scale the diagonals of the diffusion operator for the
%implicit Back-Eul operator
LImpUDiag = -dt*LSolUDiag;
LImpLDiag = -dt*LSolLDiag;
LImpMDiag = 1 - dt*LSolMDiag;



%Final variable coefficient implicit operator
Lfut = spdiags([LImpLDiag',0;LImpMDiag';0,LImpUDiag']',-1:1,Ncell,Ncell);

end