function out = dnzero(delta,alpha,Delta)
%this function returns 
%(cbar*C_M(alpha,delta)+Delta*C_Z(alpha,0,delta))^2

cbar = 11.32;

[Cx,Cm] = cxcm(alpha,delta);
Cz      = cz(alpha,0,delta);

out     =(cbar*Cm + Delta*Cz)^2;

