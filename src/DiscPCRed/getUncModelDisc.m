clc; clear;
load('trimmedF16');

%% clearing data related to Vtrim = 300
f16ss(1) = [];
Vtrims(1) = [];
utrim(:,1) = [];
xtrim(:,1) = [];

%%
Ts = 0.01;
for i=1:numel(Vtrims)
    f16ss(i).sysD = c2d(f16ss(i).sys,Ts);
end
minV = min(Vtrims);
maxV = max(Vtrims);
xV = linspace(minV,maxV,100);
% Write V in terms of z in [-1 1]
syms z1 real
V = (maxV+minV)/2 + z1*(maxV-minV)/2;
Z = (Vtrims - (maxV+minV)/2)/((maxV-minV)/2);
fV = inline(V);

%% Approximate A matrix
% =====================
nOrd = 5;
Phi = getBasis('Legendre',1,nOrd,'z');
phi = inline(Phi');
PP = zeros(numel(Vtrims),nOrd+1);
for i=1:numel(Z)
    PP(i,:) = phi(Z(i));
end

%%

for i1=1:4
    for i2 = 1:4
        for i=1:numel(Vtrims)
            a(i) = f16ss(i).sysD.A(i1,i2);
        end
        
        if var(a)<1E-6
            a = repmat(mean(a),1,numel(Vtrims));
        end
        a(abs(a)<1E-6) = 0;
        
        % Find L2 fit
        cvx_solver mosek
        cvx_precision high
        cvx_begin quiet
        variable x(nOrd+1,1)
        minimize norm(PP*x - a')
        cvx_end
        fprintf(1,'%+.4f  ',x');
        fprintf(1,'\n');
        alp(i1,i2,:) = x;
        axis tight; grid on;
    end
end

%% Approximate B matrix
% =====================
for i1=1:4
    for i2 = 1:2
        for i=1:numel(Vtrims)
            b(i) = f16ss(i).sysD.B(i1,i2);
        end
        if var(b)<1E-6
            b = repmat(mean(b),1,numel(Vtrims));
        end
        b(abs(b)<1E-6) = 0;
        subplot(4,2,(i1-1)*2+i2);
        
        % Find L2 fit
        cvx_solver mosek
        cvx_precision high
        cvx_begin quiet
        variable x(nOrd+1,1)
        minimize norm(PP*x - b')
        cvx_end
        fprintf(1,'%.5f\t',x');
        fprintf(1,'\n');
        beta(i1,i2,:) = x;
    end
end

save uncF16Disc alp beta Vtrims

%% Construct PC Model
syms z1 real

minV = min(Vtrims);
maxV = max(Vtrims);

nTerms = size(alp,3);
ord = nTerms;
Phi = getBasis('Legendre',1,nTerms-1,'z');

ns = size(beta,1);
nu = size(beta,2);

A = zeros(ns,ns);
B = zeros(ns,nu);
small = 1e-4;
for i=1:nTerms
    X = squeeze(alp(:,:,i));
    A = A + X*Phi(i);
end

for i=1:nTerms % Select terms in  B
    X = squeeze(beta(:,:,i));
    B = B + X*Phi(i);
end

%% Clean up coefficients
figure(1); clf; set(gcf,'color',[1,1,1]);
xZ = linspace(-1,1,50);
xV = linspace(minV,maxV,50);
for is=1:4
    for js = 1:4
        
        for i=1:numel(Vtrims)
            a(i) = f16ss(i).sysD.A(is,js);
        end
        if var(a)<1E-6
            a = repmat(mean(a),1,numel(Vtrims));
        end
        a(abs(a)<1E-6) = 0;
        
        figure(1);
        subplot(4,6,(is-1)*6+js); hold on;
        y = double(subs(A(is,js),xZ));
        plot(xV,y,'r',Vtrims,a,'ko','Markersize',3,'MarkerFaceColor','w'); grid on;
        str = sprintf('$A_{%d%d}$',is,js);
        ylabel(str,'Interpreter','Latex','FontSize',14);
        xlabel('$V$ (ft/s)','Interpreter','Latex','FontSize',14);
        ax = gca;
        ax.YAxis.TickLabelFormat = '%,.2f';
        axis tight;
    end
    
    for js = 1:2

        subplot(4,6,(is-1)*6+js+4); hold on;
        y = double(subs(B(is,js),xZ));
        for i=1:numel(Vtrims)
            b(i) = f16ss(i).sysD.B(is,js);
        end
        if var(b)<1E-6
            b = repmat(mean(b),1,numel(Vtrims));
            b(abs(b)<1E-6) = 0;
        end
        
        plot(xV,y,'r',Vtrims,b,'ko','Markersize',3,'MarkerFaceColor','w'); grid on;
        str = sprintf('$B_{%d%d}$',is,js);
        ylabel(str,'Interpreter','Latex','FontSize',14);
        xlabel('$V$ (ft/s)','Interpreter','Latex','FontSize',14);
        ax = gca;
        ax.YAxis.TickLabelFormat = '%,.2f';
        axis tight;
   end
end

figure(1); fig = gcf; fig.PaperPositionMode = 'auto'; 
% print('ABmatrix','-depsc','-r900');

%% Cost weight
% V(ft/s)   alp (rad)   th q  
C = [1 0 0 0; 0 1 0 0; 0 -1 1 0;]; % V alp gam
Q = C'*diag([1E-1, 1E1, 1E1])*C;

% thrust dele (deg)  
R = diag([1E-4 1E-1]);

% Q = diag([1E-2; % V (ft/s)
%     1E1;  % alp (rad)
%     1E1;  % th (rad)
%     1E-4; % q (rad/s)
%     ]);
% 
% %        thrust dele
% R = diag([1E-4   1E-1]);

save pcF16LQRModelDisc400 Q R A B



