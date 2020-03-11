
function simulateNLF16
load trimmedF16;
addpath('./F16 Model');
trimPoint = 5;
xTrim = xtrim(:,trimPoint);
uTrim = utrim(:,trimPoint);

load pcSimulationUB % Load controller
load pcF16LQRModel 

d2r = pi/180;
xp0 = [50; 0*d2r;0*d2r; 0*d2r]; % V alp th q
x0 = xp0+xTrim;
tspan = [0 10];

NMAX = numel(saveVar);
for i=1:NMAX
    [t,y] = ode45(@modifiedDynamics,tspan,[x0;0],[],saveVar(i).K,Q,R,xTrim,uTrim);
    cost(i,1) = y(end,end);
    disp(y(end,end));
end

end


function xdot = modifiedDynamics(t,x,K,Q,R,xTrim,uTrim)
xx = x(1:end-1);  % F16 states.
xp = xx - xTrim;  % Perturbation states
u = K*xp + uTrim; % total Control = Kx + uTrim.

xdotF16 = f16longi(xx,u);
xdotCost = xp'*(Q+K'*R*K)*xp; % Integral cost is based on perturbations.
xdot = [xdotF16;xdotCost];
end



