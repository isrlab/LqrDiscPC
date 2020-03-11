function J = trimCost(P)
xdot = f16longi(P(1:4),P(5:6));
Q = diag([1/100 1 1 1]);
J = xdot'*Q*xdot;
return