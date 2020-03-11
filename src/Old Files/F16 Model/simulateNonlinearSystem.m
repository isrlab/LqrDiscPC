% Simulate nonlinear dynamics with the linear controller

eqF16 = load('trimmedF16'); Vtrims = eqF16.Vtrims(1:end);
controllerInf = load('../infSynthesis_uncB1');
controllerSurr = load('../surrSynthesis_uncB1');


controller = controllerInf; col = 'r';
nTrims = numel(Vtrims);
d2r = pi/180;
dx0 = [0; 0;10*d2r; 0];
figure(1); clf;
figure(2); clf;
K = controller.saveVar(controller.NMAX).K;
yLab = {'$V$ (ft/s)','$\alpha$ (rad)','$\gamma$ (rad)'};
uLab = {'$T$ (lb)','$\delta_e$ (deg)'};

minV = min(Vtrims);
maxV = max(Vtrims);
del = (Vtrims - (maxV+minV)/2)/((maxV-minV)/2);
nu = 2;

for i=1:nTrims
    xTrim = eqF16.xtrim(:,i);
    uTrim = eqF16.utrim(:,i);
    x0 = xTrim + dx0;
    [t,x] = ode45(@f16LongiDynamics,[0 20],x0,[],xTrim,uTrim,K);
    gam = x(:,3)-x(:,2);
    y = [x(:,1) x(:,2) gam]; ny = size(y,2);
    yTrim = [xTrim(1) xTrim(2) (xTrim(3)-xTrim(2))];
    for k=1:ny
        figure(1);
        subplot(ny,1,k); hold on;
        h=plot(t,y(:,k)-yTrim(k),col,'Linewidth',1); grid on;
        %scale = 0.5*del(i)+0.5;
        %h.Color = h.Color*scale;
        
        ylabel(yLab{k},'Interpreter','Latex','Fontsize',12);
        xlabel('Time (s)','Interpreter','Latex','Fontsize',12);
    end
    
    % Controller
    nt = numel(t);
    xxTrim = repmat(xTrim',nt,1);
    uuTrim = repmat(uTrim',nt,1);
    dx = x-xxTrim;
    du = dx*K';
    u = uuTrim + du;
    for k=1:nu
        figure(2);
        subplot(nu,1,k); hold on;
        h=plot(t,u(:,k),col,'Linewidth',1); grid on;
        ylabel(uLab{k},'Interpreter','Latex','Fontsize',12);
        xlabel('Time (s)','Interpreter','Latex','Fontsize',12);
    end
end
figure(1); print -depsc ../nlSimInfState.eps
figure(2); print -depsc ../nlSimInfControl.eps


