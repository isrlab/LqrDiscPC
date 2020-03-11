% Code for synthesis Reduced Order Controller
clc; 
clear;
% Use file pcF16LQRModelDisc400 for data from system
% with Vtrims 400:100:900;
load pcF16LQRModelDisc400
z1 = sym('z1');
assume(z1,'real');

[ns,nu] = size(B);
PDF = 0.5;
alpha = 1;
filename = 'redOrdSyn';
%% Polynomial Chaos Approach
NMAX = 7; % Expansion order.
for N=1:NMAX
% for N = NMAX
    syms z1 real
    Phi = getBasis('Legendre',1,N,'z');
    Phin = kron(Phi,eye(ns));
    Phim = kron(Phi,eye(nu));

    nBasis = numel(Phi);
    NN = ns*(nBasis);
    
    M1 = innerProduct(PDF*Phi*Phi',z1,[-1,1]);
    A1 = innerProduct(PDF*kron(Phi*Phi',A),z1,[-1,1]);
    B1 = innerProduct(PDF*kron(Phi*Phi',B),z1,[-1,1]);
    
    Apc = (kron(double(M1),eye(ns,ns)))\double(A1); 
    Bpc = (kron(double(M1),eye(ns,ns)))\double(B1)*alpha;
        
    Qpc = double(innerProduct(PDF*Phin*Q*Phin',z1,[-1,1]));
    Rpc = double(innerProduct(PDF*Phim*R*Phim',z1,[-1,1]));
    
    %% Controller Synthesis
    % CVX Code
    tic;
%     cvx_solver mosek
    cvx_begin sdp quiet
    variable Ppc(NN,NN) symmetric
    
    maximize trace(Ppc)
    [(Apc'*Ppc*Apc + Qpc - Ppc) (Bpc'*Ppc*Apc)'; 
     (Bpc'*Ppc*Apc) (Rpc + Bpc'*Ppc*Bpc) ] >= 0
    cvx_end
    disp(cvx_status);
    saveVar(N).Tpc = toc;
    saveVar(N).Ppc = Ppc;
        
    % Recover controller
    Z = (Rpc + Bpc'*Ppc*Bpc)\Bpc'*Ppc*Apc;
    cvx_begin sdp quiet
    variable K(nu,ns) 
    variable X(NN,NN) symmetric
    
    Kpc = kron(eye(nBasis),K);
    invR = inv(Rpc + Bpc'*Ppc*Bpc);
    invR = (invR + invR')/2; % removing round off error for symmetry
    
    [X (Kpc + Z)';
     (Kpc + Z)  invR] >= 0
    minimize trace(X)
    cvx_end
    disp(cvx_status);
    saveVar(N).K = K;
    
    fprintf(1,'PC Order: %d, CPU Time: %.3f, trace(Ppc): %.3f, norm(K):%.3f\n',N,saveVar(N).Tpc ,trace(Ppc),norm(K));    

end
str = string(strcat('save',{' '},filename,' saveVar NMAX Apc Bpc'));
eval(str);

