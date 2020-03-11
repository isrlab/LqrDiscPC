function xdot = f16longi(x,u)
% This is  the 4-state model (V,alpha,theta,q).
% u = T,el
%
% Functions called:
% atmos()  to get the density and mach number.
% damping(),cxcm(),clcn(),cz()  do the table lookups.
%
% no engine dynamics, just T as input

xdot = zeros(4,1);

% constants
% =========
g    = 32.17;          % gravity, ft/s^2
m    = 636.94;         % mass, slugs
S    = 300;            % plan area, ft^2
cbar = 11.32;          % mean aero chord, ft
xcgr = 0.35*cbar;      % reference center of gravity
xcg  = 0.30*cbar;      % center of gravity.
h    = 10000;          % Height in ft 
r2d  = 180/pi;         % radians to degrees

Jy  = 55814;           % 1/c7;        %slug-ft^2     
Jxz = 982;             % Jy*c6;       %slug-ft^2      
Jz  = 63100;           % c3*Jxz/c4;   %slug-ft^2
Jx  = 9496;            % c9*Jxz/c4;   %slug-ft^2

% States
% =======
vt    = x(1);
alpha = x(2)*r2d; % in degrees
theta = x(3);     % Pitch angle
Q     = x(4);     % Pitch Rate --- pitching moment is M

sintheta = sin(theta);
costheta = cos(theta);
tantheta = tan(theta);

ca = cos(x(2));
sa = sin(x(2));

% Control inputs
% ===============

T     = u(1);  % Thrust
el    = u(2);  % Elevator setting in degrees.

[mach,qbar] = atmos(h,vt); %sets dynamic pressure and mach number
[Cxq,Cyr,Cyp,Czq,Clr,Clp,Cmq,Cnr,Cnp] = damping(alpha);
[Cx,Cm] = cxcm(alpha,el);
[Cl,Cn] = clcn(alpha,0);
Cz = cz(alpha,0,el);

My = qbar*S*cbar*(Cm + (.5/vt)*cbar*Cmq*Q + ...
  ((xcgr-xcg)/cbar)*(Cz+cbar*(.5/vt)*Czq*Q));


% Vdot
xdot(1) = (1/m)*(ca*(-m*g*sintheta + T + qbar*S*(Cx+cbar*(.5/vt)*Cxq*Q)) + ...
           sa*(m*g*costheta + qbar*S*(Cz+cbar*(.5/vt)*Czq*Q)));

% Alpha dot
xdot(2) = (Q - ...
  (sa/(m*vt))*(-m*g*sintheta+ T +qbar*S*(Cx+cbar*(.5/vt)*Cxq*Q))+...
  (ca/(m*vt))*(m*g*costheta + qbar*S*(Cz+cbar*(.5/vt)*Czq*Q)));

% Theta dot
xdot(3) = Q;

% Q dot
xdot(4) = My/Jy;

