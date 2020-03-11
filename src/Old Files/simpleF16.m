clc; clear;
syms z1 real

beta = 0.4;
A = [0.1658 -13.1013 -7.2748*(1 + beta*z1) -32.1739 0.2780;
    0.0018 -0.1301 0.9276*(1 + beta*z1) 0 -0.0012;
    0 -0.6436 -0.4763 0 0;
    0 0 1 0 0;
    0 0 0 0 -1;];

B = [0 -0.0706;
    0  -0.0004;
    0  -0.0157;
    0   0;
    64.94 0];

% V   alp    q    th   Thrust
C = [0 -1 0 1 0; 1 0 0 0 0];
Q = 100*C'*diag([1E2, 1])*C;

%     throttle dele (deg)  
R = diag([1E2 1E-1]);

save simpleF16 A B Q R