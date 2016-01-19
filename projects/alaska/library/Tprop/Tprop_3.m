 function [K,rho,Cp,Cpave,fracLiq] = Tprop_3(T,materialC)

% Thermal properties at CV grid points.
% Tprop version '3'.
% Saturated permafrost using Met-station derived SMCC curve.

% ________________________________________________________________

 global mlib_loc

 mlib_loc = '~/thermal/material_lib';   % location of material library

% unpack material parameters

 phi  = materialC{3};
 Kg25 = materialC{4};

% define component densities

 rhog = 2650;   % density of argilleceous mineral grains
 rhow = 1000;   % density of water

% find unfrozen water content w_u

 fname = 'AK102_smcc_2009';
 load([mlib_loc '/SMCC/' fname])
 phi_U = phi_u;       % store phi_u from file
 phi_u = interp1(-theta,phi_U,T,'linear');    % interpolate onto T grid
 phi_m = max(phi_U);                          % porosity at Met Site
 w_u   = rhow/rhog * phi_u/(1-phi_m);         % unfrozen water content

% convert w_u to volume fraction of water phi_u

 phi_u    = (rhog/rhow) * (1 - phi) .* w_u;
 L        = phi_u > phi;
 phi_u(L) = phi(L);           % phi_u cannot exceed phi
 L        = T > 0;
 phi_u(L) = phi(L);           % phi_u = phi when T > 0
 fracLiq  = phi_u ./ phi;     % fraction of water in the liq state

% find bulk density, heat capacity, and thermal conductivity

 rho = rho_permafrost( T,phi,phi_u);
 C   = Csub_permafrost(T,phi,phi_u);
 K   = Ksub_permafrost(T,phi,phi_u,Kg25);

% find the specific heat

 Cp = C ./ rho;

% Still need to correctly calculate Cpave

 Cpave = Cp;
