function CZ=cz(alfa,beta,dele)
% CZ=cz(alfa,beta,dele)


A=[.770 .241 -.100 -.416 -.731 -1.053 -1.366 -1.646 -1.917 -2.120 -2.248 -2.229];
s=.2*alfa;
k=fix(s);
if k<=-2
  k=-1;
elseif k>=9
  k=8;
end
da=s-k;
L=k+fix(1.1*sign(da));
k=k+3;
L=L+3;
s=A(k)+abs(da)*(A(L)-A(k));
CZ=s*(1-(beta/57.3)^2)-.19*dele/25;

