% Simulate nonlinear dynamics with the linear controller
clear;
eqF16 = load('/home/vish0908/Documents/Github/ACC_LQRDiscPC/src/Old Files/F16 Model/trimmedF16'); 
%% clearing data related to Vtrim = 300
% eqF16.f16ss(1) = [];
% eqF16.Vtrims(1) = [];
% eqF16.utrim(:,1) = [];
% eqF16.xtrim(:,1) = [];
%%
Vtrims = eqF16.Vtrims(1:end);
% controllerInf = load('../infSynthesis_uncB1');
controllerSurr = load('redOrd300Syn1_3.mat');


controller = controllerSurr; col = 'b';
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
    [t,x] = ode45(@f16LongiDynamics,[0 1800],x0,[],xTrim,uTrim,K);
    gam = x(:,3)-x(:,2);
    y = [x(:,1) x(:,2) gam]; ny = size(y,2);
    yTrim = [xTrim(1) xTrim(2) (xTrim(3)-xTrim(2))];
    for k=1:ny
        figure(1);
        subplot(ny,1,k); hold on;
        h=plot(t,y(:,k)-yTrim(k),col,'Linewidth',1); grid on;
        %scale = 0.5*del(i)+0.5;
        %h.Color = h.Color*scale;
        titleText = string(strcat('minV = ',{' '},num2str(minV),' ft/s'));
        title(titleText);
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
        titleText = string(strcat('minV = ',{' '},num2str(minV),' ft/s'));
        title(titleText);
        ylabel(uLab{k},'Interpreter','Latex','Fontsize',12);
        xlabel('Time (s)','Interpreter','Latex','Fontsize',12);
    end
end
% figure(1); print -depsc nlSimRed13State300.eps
% figure(2); print -depsc nlSimRed13Control300.eps

% figure(1); print -djpeg nlSimRed13State300
% figure(2); print -djpeg nlSimRed13Control300
