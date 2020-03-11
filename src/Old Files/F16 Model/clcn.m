function [CL,CN]=clcn(alfa,beta)
%[CL,CN]=clcn(alfa,beta)

AL=[0     0     0     0     0     0     0     0     0     0     0     0
-.001 -.004 -.008 -.012 -.016 -.019 -.020 -.020 -.015 -.008 -.013 -.015
-.003 -.009 -.017 -.024 -.030 -.034 -.040 -.037 -.016 -.002 -.010 -.019
-.001 -.010 -.020 -.030 -.039 -.044 -.050 -.049 -.023 -.006 -.014 -.027
    0 -.010 -.022 -.034 -.047 -.046 -.059 -.061 -.033 -.036 -.035 -.035
 .007 -.010 -.023 -.034 -.049 -.046 -.068 -.071 -.060 -.058 -.062 -.059
 .009 -.011 -.023 -.037 -.050 -.047 -.074 -.079 -.091 -.076 -.077 -.076];
%
AN=[0    0    0    0    0    0    0    0    0     0     0     0
 .018 .019 .018 .019 .019 .018 .013 .007 .004 -.014 -.017 -.033
 .038 .042 .042 .042 .043 .039 .030 .017 .004 -.035 -.047 -.057
 .056 .057 .059 .058 .058 .053 .032 .012 .002 -.046 -.071 -.073
 .064 .077 .076 .074 .073 .057 .029 .007 .012 -.034 -.065 -.041
 .074 .086 .093 .089 .080 .062 .049 .022 .028 -.012 -.002 -.013
 .079 .090 .106 .106 .096 .080 .068 .030 .064  .015  .011 -.001];
%
s=.2*alfa;
k=fix(s);
if k<=-2
  k=-1;
elseif k>=9
  k=8;
end
da=s-k;
L=k+fix(1.1*sign(da));
s=.2*abs(beta);
m=fix(s);
if m==0
  m=1;
elseif m>=6
  m=5;
end
db=s-m;
n=m+fix(1.1*sign(db));
k=k+3;
L=L+3;
m=m+1;
n=n+1;
%
% CL
%
t=AL(m,k);
u=AL(n,k);
v=t+abs(da)*(AL(m,L)-t);
w=u+abs(da)*(AL(n,L)-u);
dum=v+(w-v)*abs(db);
CL=dum*sign(beta);
%
% CN
%
t=AN(m,k);
u=AN(n,k);
v=t+abs(da)*(AN(m,L)-t);
w=u+abs(da)*(AN(n,L)-u);
dum=v+(w-v)*abs(db);
CN=dum*sign(beta);

