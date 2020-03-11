clear;clc;
syms z1 real
N = 6;
pdf = 0.5;
Phi = getBasis('Legendre',1,N,'z');

nBasis = numel(Phi);
for i=1:nBasis
    for j=1:nBasis
        E(i,j) = innerProduct(pdf*Phi(i)^2*Phi(j),z1,[-1,1]);
    end
end


