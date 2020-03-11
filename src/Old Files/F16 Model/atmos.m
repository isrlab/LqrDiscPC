function [mach,qbar]=atmos(alt,vt)
%This function 
%[mach,qbar]=atmos(alt,vt)
%Calculates the mach number and dynamic pressure for a given
%altitude and total velocity.


%
% constants
%
rho0=2.377e-3;
%
% temperature
%
tfac=1-.703e-5*alt;
temp=519*tfac;
if alt>=35000
  temp=390;
end
%
% density
%
rho=rho0*tfac^4.14;
%
% Mach number
%
mach=vt/sqrt(1.4*1716.3*temp);
%
% pressures
%
qbar=.5*rho*vt*vt;
% ps=1715.0*rho*temp;
