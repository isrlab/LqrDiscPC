clc; clear all;

load pcF16LQRModel; uncB = 1;
%load simpleF16; uncB = 0;
ns = size(A,1);
In = eye(ns);

z1 = sym('z1');
assume(z1,'real');

[ns,nu] = size(B);
PDF = 0.5;

%% Polynomial Chaos Approach
NMAX = 8;
for N=1:NMAX
    Phi = getBasis('Legendre',1,N,'z');
    Phin = kron(Phi,eye(ns));
    Phim = kron(Phi,eye(nu));
    
    nBasis = numel(Phi);
    NN = ns*(nBasis);
 
    count = 1; th = sym([]);
    for i=1:nBasis
        for j=i:nBasis
            if i==j
                th(count,1) = Phi(i)^2;
            else
                th(count,1) = 2*Phi(i)*Phi(j);
            end
            count = count + 1;
        end
    end
    thn = kron(th,In);
    
    Qbar = double(innerProduct(PDF*Phin*Q*Phin',z1,[-1,1]));
   
    % CPU Time experiments
    % ====================
    %     tic;
    %     [C,IA,IC] = unique(e1);
    %     tt1 = double(innerProduct(PDF*C,z1,[-1,1]));
    %     E1 = reshape(tt1(IC),size(e1));
    %     toc;

    % First inner product
    e1 = kron(Phi*Phi',A'*thn');
    E1 = double(innerProduct(PDF*e1,z1,[-1,1]));

    % Second inner-product
    e2 = kron(kron(Phi,th),B); 
    
    h = waitbar(0,'Please wait...');
    XX = PDF*e2*inv(R)*e2.'; % Square matrix.
    nx = size(XX,1); 
    count = 1;
    countMax = nx*(nx+1)/2; Z1 = zeros(nx,nx);
    for i=1:nx
        tmpV1 = [];
        parfor j=i:nx
            tmpV1(j,1) = double(innerProduct(XX(i,j),z1,[-1,1]));
        end
        
        tmpV2 = tmpV1(i:nx);
        
        Z1(i,i:nx) = tmpV2';
        Z1(i:nx,i) = tmpV2;
        waitbar(i/nx,h,sprintf('N: %d, Completed: %i',N,i));
    end
    close(h);
    
    Z = sqrtm(Z1);
    W =  double(innerProduct(PDF*Phi.*Phi,z1,[-1,1]));
    M = nBasis*(nBasis+1)/2;
    disp('Starting cvx ...');
    % CVX Code
    cvx_clear;
    
    tic;
    cvx_solver sedumi
    %cvx_precision best
    cvx_begin sdp quiet
    
    variable bigP(NN,NN) symmetric
    expression costP(ns,ns)
    
    bigP >= 0
        
    % Construct Ppc from bigP
    Ppc = [];
    for i=1:nBasis
        i1 = (i-1)*ns+1; i2 = i*ns;
        costP = costP + W(i)*bigP(i1:i2,i1:i2);
        for j=i:nBasis
            j1 = (j-1)*ns+1; j2 = j*ns;
            Ppc = [Ppc;bigP(i1:i2,j1:j2)]; % Not the best way ... but ...
        end
    end

    PPpc = kron(eye(nBasis),Ppc);

    [(E1*PPpc + PPpc'*E1' + Qbar) (PPpc'*Z1);
       Z1*PPpc      eye(ns*M*nBasis)] >=0
    maximize trace(costP)
    cvx_end
    Tpc = toc;
    fprintf(1,'cvx status: %s, cost: %f, cpuTime: %f ',cvx_status,trace(costP),Tpc);
    saveVar(N).P = bigP;
    saveVar(N).Tpc = Tpc;
    saveVar(N).cost = trace(costP);
    
    %% Recover controller
    %disp('Controller recovery ...');
    Pp = Phin'*bigP*Phin;
    Rbar = double(innerProduct(PDF*Phim*R*Phim',z1,[-1,1]));
    M1 = innerProduct(PDF*Phin*Pp*B*Phim',z1,[-1,1]);
    M2 = innerProduct(PDF*Phin*Pp*B*inv(R)*B'*Pp*Phin',z1,[-1,1]);
    
    M1 = double(M1);
    M2 = double(M2);
    
    % Recover controller
    tic;
    cvx_solver mosek
    cvx_begin sdp quiet
    variable K(nu,ns)
    variable X(ns*nBasis,ns*nBasis) symmetric
    
    K1 = kron(eye(nBasis),K);
    
    X >= 0
    [(M1*K1 + K1'*M1' + M2 - X) K1';
        K1 -inv(Rbar)] <= 0
    minimize trace(X)
    cvx_end
    saveVar(N).Tpc = saveVar(N).Tpc + toc;
    saveVar(N).K = K;
    fprintf(1,'normK: %f\n',norm(K));
    disp(cvx_status);
end

str = sprintf('save infSynthesis_uncB%d saveVar NMAX',uncB);
eval(str)


