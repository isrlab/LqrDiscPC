% Plots 1a, 1b, and 2a from the paper.

clc; 
clear;
load pcF16LQRModelDisc400
Ts = 0.01;

Vtrims = 400:100:900;
nV = numel(Vtrims);
Ad = inline(A);
Bd = inline(B);

controllerSurr = load('redOrdSyn'); 

%% Closed Loop Poles
del = linspace(-1,1,nV);
isStable = zeros(nV,1);
for i=1:nV
    sys = ss(Ad(del(i)),Bd(del(i)),eye(size(A)),0,Ts);
    E(i,:) = eig(Ad(del(i)));
    max(abs(E(i,:)))
    if max(abs(E(i,:))) < 1
        isStable(i,1) = 1;
        Wc = gram(sys,'c');
        invWc = inv(Wc);
        Lw(i,1) = det(invWc);
    else
        Lw(i,1) = nan;
    end
    
    clpA = Ad(del(i)) + Bd(del(i))*controllerSurr.saveVar(end).K;
    Esur(i,:) = eig(clpA);
    
end

%% Develop surrogate model
Apc = controllerSurr.Apc;
Bpc = controllerSurr.Bpc;

%% Vtrims vs Ctrb
figure(4); clf;

semilogy(Vtrims,Lw,'ko-','Linewidth',1,'MarkerSize',5,'MarkerFaceColor','w');
ylabel('\textbf{det}$(W^{-1}_c)$','Interpreter','Latex','Fontsize',14);
xlabel('$V_\textrm{trim}$','Interpreter','Latex','Fontsize',14);
grid on
figure(4); fig = gcf; fig.PaperPositionMode = 'auto'; 
% print('ctrbEnergy','-depsc','-r900');
%% Closed Loop Poles
figure(5); clf;

Rsur = real(vec(Esur));
Isur = imag(vec(Esur));

plot(Rsur,Isur,'bo','MarkerFaceColor','b','MarkerSize',5); hold on;
viscircles([0 0],1);
grid on;
figure(5); fig = gcf; fig.PaperPositionMode = 'auto'; 
% print('clLoopPoles','-depsc','-r900');
%% Open Loop Poles
figure(6); clf;
R = real(vec(E));
I = imag(vec(E));
% 
Epc = eig(Apc);
Rpc = real(vec(Epc));
Ipc = imag(vec(Epc));

plot(R,I,'ro','MarkerFaceColor','r','MarkerSize',5); hold on;
plot(Rpc,Ipc,'bp','MarkerSize',5,'MarkerFaceColor','b');
% axis equal

grid on;
h = legend('Poles from uncertain system','Poles from reduced order model');

figure(6); fig = gcf; fig.PaperPositionMode = 'auto'; 


