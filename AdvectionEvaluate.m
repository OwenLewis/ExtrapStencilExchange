%This function evaluates the (explicit) advection terms of a
%reaction-advection-diffusion system. The scheme used is a very simple
%centered differencing scheme. Obviously, this should not be used as part
%of a Forward-Euler time stepping scheme. It will be unconditionally
%unstable. However, as part of an IMEX time integration scheme, it will be
%quite useful. 
%
% function syntax:
%
%     advectionTerms = AdvectionEvaluate(conc,veloc)
%
%
%     inputs:
%         conc is a length Ncell array which contains the concentrations
%           (at cell centeres) of the species which is being advected. 
%
%         veloc is a length Nedges+2 array which contains the velocites (at
%           cell edges) with which the species is being advected. Again,
%           the first and last entries (ghost cells) should have aleady
%           been populated with whatever values respect the B.C's
%     output:
%         advectionTerms is a length Ncell array defined at cell centers
%           which contains the explicit evaluation
%           of the advection of the species 'conc' at speed 'veloc'

function advectionTerms = AdvectionEvaluate(conc,veloc)


%Lets 'import' the two big global structs
global GelState GelSimParams


%Here are some parameters we need to define the sizes of things
hx = GelSimParams.hx;



%We should check to make sure no downwinding is occuring.
down = sum(veloc < 0);
if down > 0
    disp('Woops, you are downwinding')
end

%We need to evaluate the concentration in the "ghost cell" which
%does not exist. We can do this via linear extrapolation from the
%two left-most interior cells. 
ghostvalue = 2*conc(1) - conc(2);
conc = [ghostvalue;conc];
%Flux through this edge is velocity times the concentration in the
%cell center TO THE LEFT
flux = conc.*veloc;
%Do a centered difference of flux to get the term you want
advectionTerms = diff(flux)/hx;


end
