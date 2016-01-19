 function [K,rho,Cp,C,diagC] = Tprop_4(T,materialC,spinup)

% Thermal properties at CV grid points.
% Tprop version '4'.
% Saturated permafrost using Anderson et al (1973) unfrozen water equations.
% This version will work for coarse or fine-grained permafrost.

% For now, it is assumed the type of material is the same at all depths,
% although the porosity can vary.  Thus, imat and matfrac are the same at
% all depths.

% During the spinup phase, the pores are assumed to be filled with ice and
% unfrozen water.  After spinup, some air space may develop as ice converts
% to unfrozen water.  The total amount of water in the pore spaces is assumed
% to remain fixed with time.

%   T           temperature (C)
%   materialC   parameters specifying the type of material (cell array)
% ________________________________________________________________

 global mlib_loc

 mlib_loc = '~/thermal/material_lib';   % location of material library

% unpack cell array

 imat    = materialC{1}(1); 
 matfrac = materialC{2}(1);
 Kg25    = materialC{3};
 phi     = materialC{4};
 Mw      = materialC{5};        % values at end of spinup

% define component densities

 rhog = 2650;       % density of argilleceous mineral grains
 rhoi =  917;       % density of ice
 rhow = 1000;       % density of water

% find unfrozen water content w_u

 switch imat
 case 0             % coarse materials
   w_u    = zeros(size(T));     % assume w_u = 0 for T < 0 for coarse materials
   L      = T >= 0;
   w_u(L) = 1;
 otherwise          % fine materials
   [~,w_u] = unfrozenWater(T,imat,matfrac);
 end

% convert w_u to volume fraction of water phi_u

 phi_u    = (rhog/rhow) * (1 - phi) .* w_u;
 L        = phi_u > phi;
 phi_u(L) = phi(L);           % phi_u cannot exceed phi
 L        = T > 0;
 phi_u(L) = phi(L);           % phi_u = phi when T > 0

% find volume fractions of ice and air

 if spinup
   phi_i    = phi - phi_u;
   phi_a    = zeros(size(phi));
 else
   phi_i    = (Mw - rhow*phi_u) / rhoi;
   L        = phi_i < 0;
   phi_i(L) = 0;
   phi_a    = phi - (phi_i + phi_u);
 end

% find total mass of pore water (frozen and unfrozen)

 Mw = rhoi*phi_i + rhow*phi_u;

% find bulk density, heat capacity, and thermal conductivity

 rho = rho_permafrost(phi,phi_i,phi_u);
 C   = Csub_permafrost(T,phi,phi_i,phi_u);  % kinetic-energy component of heat capacity
 K   = Ksub_permafrost(T,phi,phi_i,phi_u,Kg25);
 Cp  = C ./ rho;

% find dphi_u/dT at T

 switch imat
 case 0
   dphiudT  = zeros(size(T));
 otherwise
   dT       = 0.05;
   [~,w_up] = unfrozenWater(T+dT, imat,matfrac);
   [~,w_um] = unfrozenWater(T-dT, imat,matfrac);

   phi_up    = (rhog/rhow) * (1 - phi) .* w_up;
   L         = phi_up > phi;
   phi_up(L) = phi(L);
   L         = T > 0;
   phi_up(L) = phi(L);

   phi_um    = (rhog/rhow) * (1 - phi) .* w_um;
   L         = phi_um > phi;
   phi_um(L) = phi(L);
   L         = T > 0;
   phi_um(L) = phi(L);

   dphiudT   = (phi_up - phi_um) / (2*dT);
 end

% store volume fractions in diagnostics cell array

 diagC = {Mw,phi_i,phi_u,phi_a,dphiudT};
