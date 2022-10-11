function xdot = fun_xdot(x,u)

global m1 m2 L1 L2 g val Nx Nu pert

xdot = [x(5:8); fun_qddot(x,u)];
%xdot = fun_qddotssp11(x,u,dt);
