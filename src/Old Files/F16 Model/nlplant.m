function xdot = nlplant(x,u)
% This is  the 6-state model (V,alpha,theta,q).
% u = T,el
%
% Functions called:
% atmos()  to get the density and mach number.
% damping(),cxcm(),clcn(),cz()  do the table lookups.
%
% no engine dynamics, just T as input

xdot = zeros(6,1);

% constants
% =========
g    = 32.17;          %gravity, ft/s^2
m    = 636.94;         %mass, slugs
S    = 300;            %plan area, ft^2
cbar = 11.32;          %mean aero chord, ft
xcgr = 0.35*cbar;      %reference center of gravity
xcg  = 0.30*cbar;      %center of gravity.
%
r2d  = 180/pi;            %radians to degrees
%
%     %NasaData        
Jy  = 55814;           %1/c7;        %slug-ft^2     
Jxz = 982;             %Jy*c6;       %slug-ft^2      
Jz  = 63100;           %c3*Jxz/c4;   %slug-ft^2
Jx  = 9496;            %c9*Jxz/c4;   %slug-ft^2

% States
% =======

npos  = x(1); % north position 
alt   = x(2); % altitude
theta = x(3); % Pitch angle
vt    = x(4);
alpha = x(5)*r2d; % in degrees
Q     = x(6);    % Pitch Rate--- pitching moment is M
sa    = sin(x(5)); %sin(alpha) 
ca    = cos(x(5)); %cos(alpha) 

% Control inputs
% ===============
T     = u(1);  % Scaled Thrust
el    = u(2);  % Elevator setting in degrees.

[mach,qbar] = atmos(alt,vt); %sets dynamic pressure and mach number


U = vt*ca;  %directional velocities.
W = vt*sa;

sintheta = sin(theta);
costheta = cos(theta);
tantheta = tan(theta);

xdot(1) = U*(costheta) + W*(sintheta);
xdot(2) = U*(sintheta) - W*(costheta);
xdot(3) = Q;


[Cxq,Cyr,Cyp,Czq,Clr,Clp,Cmq,Cnr,Cnp] = damping(alpha);
[Cx,Cm]                               = cxcm(alpha,el);
[Cl,Cn]                               = clcn(alpha,0);
Cz                                    = cz(alpha,0,el);

xdot(4) = (1/m)*(...
  ca*(-m*g*sintheta + T + qbar*S*(Cx+cbar*(.5/vt)*Cxq*Q)) + ...
  sa*(m*g*costheta + qbar*S*(Cz+cbar*(.5/vt)*Czq*Q)));


xdot(5) = (Q - ...
  (sa/(m*vt))*(-m*g*sintheta+ T +qbar*S*(Cx+cbar*(.5/vt)*Cxq*Q))+...
  (ca/(m*vt))*(m*g*costheta + qbar*S*(Cz+cbar*(.5/vt)*Czq*Q)));
out(5) = alphadot;

My = qbar*S*cbar*(Cm + (.5/vt)*cbar*Cmq*Q + ...
  ((xcgr-xcg)/cbar)*(Cz+cbar*(.5/vt)*Czq*Q));

xdot(6) = My/Jy;


