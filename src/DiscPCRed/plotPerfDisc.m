function plotPerfDisc
% Function to display plots 2b and 3 from the paper.

controllerSurr = load('redOrdSyn'); 
plant = load('pcF16LQRModelDisc400');

%% Controller gain vs approximation order
figure(7); clf; hold on;

plot(1:controllerSurr.NMAX,getControllerNorm(controllerSurr),'k','Linewidth',1);
xlabel('PC Approximation Order','Fontsize',12,'Interpreter','Latex');
ylabel('$\|\mathbf{K}\|_2$','Fontsize',12,'Interpreter','Latex');
grid on; box on;
figure(7); fig = gcf; fig.PaperPositionMode = 'auto'; 
% print('kGain','-depsc','-r900');

%% Plot system trajectories.

d2r = pi/180;
x0 = [0; 0; 30*d2r; 0];
xLab = {'$V$ (ft/s)','$\alpha$ (rad)','$\theta$ (rad)','$q$ (rad/s)'};

figure(10); clf;
orderAppx = controllerSurr.NMAX; t = 50/0.01;
plotSystemTraj(controllerSurr,plant,orderAppx,t,'b',x0);
figure(10); fig = gcf; fig.PaperPositionMode = 'auto'; 
% print('yTrajRed','-depsc','-r900');
end

function normK = getControllerNorm(controller)
for i=1:controller.NMAXg
    normK(i,1) = norm(controller.saveVar(i).K);
end
end


function plotSystemTraj(controllerInf,plant,orderAppx,t,col,x0)
K = controllerInf.saveVar(orderAppx).K;
ACLP = inline(plant.A + plant.B*K);
yLab = {'$V$ (ft/s)','$\alpha$ (rad)','$\gamma$ (rad)'};
nSamp = 50;
del = linspace(-1,1,nSamp);
[ns,nu] = size(plant.B);
ny = 3;
for i=1:nSamp
    Aclp = ACLP(del(i));
    [t1,x] = dynamicsDisc(t,x0,Aclp);
    V = x(1,:)';
    alp = x(2,:)';
    gam = -x(2,:)' + x(3,:)';
    y = [V alp gam];
    Del = ones(size(t1))*del(i);
    for j=1:ny
        subplot(ny,1,j); hold on; h = plot(t1,y(:,j),col,'Linewidth',1); axis tight
        scale = 0.5*del(i)+0.5;
        h.Color = h.Color*scale;
    end
end
for j=1:ny
    subplot(ny,1,j);
    ylabel(yLab{j},'Interpreter','Latex','Fontsize',12);
    if(j==3)
        xlabel('Time (s)','Interpreter','Latex','Fontsize',12);
    end
    grid on;
end
end

function [tSim,xDisc] = dynamicsDisc(nT,x0,A)
    xDisc = zeros(size(x0,1),nT);
    xDisc(:,1) = x0;
    for i=1:nT-1
        xDisc(:,i+1) = A*xDisc(:,i);
    end
    tSim = 0.01*[1:nT]';
end

function [h] = plotResults(NMAX,trP,symb1,symb2)
h = semilogy(1:NMAX,trP,'k','Linewidth',1); hold on;
h = semilogy(1:NMAX,trP,symb1,'MarkerFace','w');
end

function trP = getControllerPerformanceData(controllerDataFile,plantDataFile)
controller = load(controllerDataFile);
plant = load(plantDataFile);
[ns,nu] = size(plant.B);
trP = getTrace(controller,plant.A,plant.B,plant.Q,plant.R);
end


function trP = getTrace(dataVar,A,B,Q,R)
[ns,nu ] =size(B);
for i=1:dataVar.NMAX
    Aclp = A+B*dataVar.saveVar(i).K;
    ACLP = inline(Aclp);

    Tspan = [0 300];
    
    param = linspace(-1,1,100);
    tmpVar = zeros(numel(param),1);
    
    parfor j=1:numel(param)
        [t,y] = ode45(@dynamicsP,Tspan,zeros(ns^2,1),[],ACLP(param(j)),Q,R,dataVar.saveVar(i).K);
        P = reshape(y(end,:),[ns,ns]);
        tmpVar(j) = trace(P);
    end
    trP(i) = mean(tmpVar);
end
end

function xdot = dynamicsP(t,x,A,Q,R,K)
    Sdot = (A^t)'*(Q+K'*R*K)*(A^t);
    xdot = vec(Sdot);
end







