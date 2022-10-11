function xdot = fc(x,u)

m1=2;
m2=3;
L1=2;L2=3;
r1=1;r2=1.5;
g = 9.81; % gravity




xdot = [x(5); x(6);x(7); x(8); fun_qddot(x,u)];