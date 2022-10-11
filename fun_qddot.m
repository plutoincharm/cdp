function qddot = fun_qddot(x,u) 
m1=50.172;
m2=7.4;
m3=3.411;
L2=0.4418;L3=0.4033;
r2=0.1872;r3=0.1738;
g = 9.81; % gravity


MI2 = (1/3).*m2.*L2.^2;
MI3 = (1/3).*m3.*L3.^2;
x1 = x(1);
y1 = x(2);
tht2 = x(3);
tht3 = x(4);
vx1 = x(5);
vy1 = x(6);
omg2 = x(7);
omg3 = x(8);
T1 = u(1);
T2 = u(2);
F1 = u(3);
F2 = u(4);



Mmat=[m2+m3,0,L2.*m3.*sin(tht2)+m2.*r2.*sin(tht2),m3.*r3.*sin(tht3);0, ...
  m2+m3,(-1).*L2.*m3.*cos(tht2)+(-1).*m2.*r2.*cos(tht2),(-1).*m3.* ...
  r3.*cos(tht3);L2.*m3.*sin(tht2)+m2.*r2.*sin(tht2),(-1).*L2.*m3.* ...
  cos(tht2)+(-1).*m2.*r2.*cos(tht2),L2.^2.*m3+MI2+m2.*r2.^2,L2.*m3.* ...
  r3.*cos(tht2+(-1).*tht3);m3.*r3.*sin(tht3),(-1).*m3.*r3.*cos(tht3) ...
  ,L2.*m3.*r3.*cos(tht2+(-1).*tht3),MI3+m3.*r3.^2];



phi=[F1+(-1).*L2.*m3.*omg2.^2.*cos(tht2)+(-1).*m2.*omg2.^2.*r2.*cos( ...
  tht2)+(-1).*m3.*omg3.^2.*r3.*cos(tht3),F2+(-1).*g.*m1+(-1).*g.*m2+ ...
  (-1).*g.*m3+(-1).*L2.*m3.*omg2.^2.*sin(tht2)+(-1).*m2.*omg2.^2.* ...
  r2.*sin(tht2)+(-1).*m3.*omg3.^2.*r3.*sin(tht3),T1+(-1).*T2+g.*L2.* ...
  m3.*cos(tht2)+g.*m2.*r2.*cos(tht2)+(-1).*L2.*m3.*omg3.^2.*r3.*sin( ...
  tht2+(-1).*tht3),T2+g.*m3.*r3.*cos(tht3)+L2.*m3.*omg2.^2.*r3.*sin( ...
  tht2+(-1).*tht3)];



qddot = inv(Mmat)*phi'

