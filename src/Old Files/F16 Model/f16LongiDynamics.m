function xdot = f16LongiDynamics(t,x,xtrim,utrim,K)
u = utrim + K*(x-xtrim);
xdot = f16longi(x,u);
%disp(xdot');

