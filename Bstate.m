function B = Bstate(x,u,dt)

global m1 m2 L1 L2  g
b11=(1/2).*dt.^2.*L2.^2.*m2.*(L1.^2.*L2.^2.*m1.*m2+L1.^2.*L2.^2.* ...
  m2.^2+(-1).*L1.^2.*L2.^2.*m2.^2.*cos(x(1)+(-1).*x(2)).^2).^(-1);

b12=(1/2).*dt.^2.*((-1).*L2.^2.*m2.*(L1.^2.*L2.^2.*m1.*m2+L1.^2.* ...
  L2.^2.*m2.^2+(-1).*L1.^2.*L2.^2.*m2.^2.*cos(x(1)+(-1).*x(2)).^2) ...
  .^(-1)+(-1).*L1.*L2.*m2.*cos(x(1)+(-1).*x(2)).*(L1.^2.*L2.^2.*m1.* ...
  m2+L1.^2.*L2.^2.*m2.^2+(-1).*L1.^2.*L2.^2.*m2.^2.*cos(x(1)+(-1).* ...
  x(2)).^2).^(-1));

b21=(-1/2).*dt.^2.*L1.*L2.*m2.*cos(x(1)+(-1).*x(2)).*(L1.^2.*L2.^2.* ...
  m1.*m2+L1.^2.*L2.^2.*m2.^2+(-1).*L1.^2.*L2.^2.*m2.^2.*cos(x(1)+( ...
  -1).*x(2)).^2).^(-1);

b22=(1/2).*dt.^2.*(L1.^2.*(m1+m2).*(L1.^2.*L2.^2.*m1.*m2+L1.^2.* ...
  L2.^2.*m2.^2+(-1).*L1.^2.*L2.^2.*m2.^2.*cos(x(1)+(-1).*x(2)).^2) ...
  .^(-1)+L1.*L2.*m2.*cos(x(1)+(-1).*x(2)).*(L1.^2.*L2.^2.*m1.*m2+ ...
  L1.^2.*L2.^2.*m2.^2+(-1).*L1.^2.*L2.^2.*m2.^2.*cos(x(1)+(-1).* ...
  x(2)).^2).^(-1));

b31=(-1/16).*dt.*L1.^(-4).*L2.^(-3).*m2.^(-1).*(m1+m2.*sin(x(1)+(-1).* ...
  x(2)).^2).^(-2).*(m1+m2.*sin((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+( ...
  -1).*x(2)).^2).^(-1).*(4.*L1.*((-4).*L1.*L2.^3.*m2.*(m1+m2.*sin( ...
  x(1)+(-1).*x(2)).^2).^2+(-4).*dt.*L1.*L2.^3.*m2.^2.*x(4).*cos( ...
  x(1)+(-1).*x(2)).*(m1+m2.*sin(x(1)+(-1).*x(2)).^2).*sin((1/2).* ...
  dt.*(x(3)+(-1).*x(4))+x(1)+(-1).*x(2))+(-2).*dt.^2.*L2.*m2.*cos( ...
  x(1)+(-1).*x(2)).*sin((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+(-1).* ...
  x(2)).*(L2.*m2.*cos(x(1)+(-1).*x(2)).*((-1).*u(1)+u(2)+L1.*(g.*m1.* ...
  sin(x(1))+m2.*(g.*sin(x(1))+L2.*x(4).^2.*sin(x(1)+(-1).*x(2)))))+ ...
  L1.*(m1+m2).*(u(2)+(-1).*L2.*m2.*((-1).*L1.*x(3).^2.*sin(x(1)+(-1).* ...
  x(2))+g.*sin(x(2))))))+m2.*cos((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+( ...
  -1).*x(2)).*(((-2).*dt.*L2.^2.*m2+(-2).*dt.*L1.*L2.*m2.*cos(x(1)+( ...
  -1).*x(2))).*(2.*dt.*L2.*((-1).*u(1)+u(2))+2.*dt.*L1.*(u(2).*cos(x(1)+( ...
  -1).*x(2))+L2.^2.*m2.*x(4).^2.*sin(x(1)+(-1).*x(2))+g.*L2.*(m1.* ...
  sin(x(1))+m2.*cos(x(2)).*sin(x(1)+(-1).*x(2))))+L1.^2.*L2.*x(3).*( ...
  (-4).*m1+m2.*((-4).*sin(x(1)+(-1).*x(2)).^2+dt.*x(3).*sin(2.*( ...
  x(1)+(-1).*x(2)))))).*sin((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+(-1).* ...
  x(2))+(-2).*dt.*L2.*(2.*dt.*L2.^2.*m2.*((-1).*u(1)+u(2))+2.*dt.* ...
  L1.^3.*L2.*m2.*(m1+m2).*x(3).^2.*sin(x(1)+(-1).*x(2))+dt.*L1.*L2.* ...
  m2.*((-2).*(u(1)+(-2).*u(2)).*cos(x(1)+(-1).*x(2))+2.*L2.^2.*m2.* ...
  x(4).^2.*sin(x(1)+(-1).*x(2))+2.*g.*L2.*(m1.*sin(x(1))+m2.*cos( ...
  x(2)).*sin(x(1)+(-1).*x(2))))+L1.^2.*(2.*dt.*(m1+m2).*u(2)+2.*dt.* ...
  g.*L2.*m2.*(m1+m2).*cos(x(1)).*sin(x(1)+(-1).*x(2))+L2.^2.*m2.*(( ...
  -4).*m1.*(x(3)+(-1).*x(4))+m2.*((-4).*x(3).*sin(x(1)+(-1).*x(2)) ...
  .^2+4.*x(4).*sin(x(1)+(-1).*x(2)).^2+dt.*x(3).^2.*sin(2.*(x(1)+( ...
  -1).*x(2)))+dt.*x(4).^2.*sin(2.*(x(1)+(-1).*x(2))))))).*sin((1/2) ...
  .*dt.*(x(3)+(-1).*x(4))+x(1)+(-1).*x(2))+(-2).*L1.*((-1).*dt.*L2.* ...
  m2.*cos(x(1)+(-1).*x(2)).*(2.*dt.*L2.*((-1).*u(1)+u(2))+2.*dt.*L1.*( ...
  u(2).*cos(x(1)+(-1).*x(2))+L2.^2.*m2.*x(4).^2.*sin(x(1)+(-1).*x(2))+ ...
  g.*L2.*(m1.*sin(x(1))+m2.*cos(x(2)).*sin(x(1)+(-1).*x(2))))+ ...
  L1.^2.*L2.*x(3).*((-4).*m1+m2.*((-4).*sin(x(1)+(-1).*x(2)).^2+dt.* ...
  x(3).*sin(2.*(x(1)+(-1).*x(2)))))).*sin((1/2).*dt.*(x(3)+(-1).* ...
  x(4))+x(1)+(-1).*x(2))+(-2).*dt.*L2.*sin((1/2).*dt.*(x(3)+(-1).* ...
  x(4))+x(1)+(-1).*x(2)).*(2.*L1.*L2.^2.*m2.*x(4).*(m1+m2.*sin(x(1)+ ...
  (-1).*x(2)).^2)+dt.*L2.*m2.*cos(x(1)+(-1).*x(2)).*((-1).*u(1)+u(2)+ ...
  L1.*(g.*m1.*sin(x(1))+m2.*(g.*sin(x(1))+L2.*x(4).^2.*sin(x(1)+(-1) ...
  .*x(2)))))+dt.*L1.*(m1+m2).*(u(2)+(-1).*L2.*m2.*((-1).*L1.*x(3).^2.* ...
  sin(x(1)+(-1).*x(2))+g.*sin(x(2))))))));


b32=(-1/16).*dt.*L1.^(-4).*L2.^(-3).*m2.^(-1).*(m1+m2.*sin(x(1)+(-1).* ...
  x(2)).^2).^(-2).*(m1+m2.*sin((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+( ...
  -1).*x(2)).^2).^(-1).*(4.*L1.*(4.*L1.*L2.^3.*m2.*(m1+m2.*sin(x(1)+ ...
  (-1).*x(2)).^2).^2+4.*dt.*L1.*L2.^2.*m2.*x(4).*(L1.*(m1+m2)+L2.* ...
  m2.*cos(x(1)+(-1).*x(2))).*(m1+m2.*sin(x(1)+(-1).*x(2)).^2).*sin(( ...
  1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+(-1).*x(2))+2.*dt.^2.*(L1.*(m1+ ...
  m2)+L2.*m2.*cos(x(1)+(-1).*x(2))).*sin((1/2).*dt.*(x(3)+(-1).* ...
  x(4))+x(1)+(-1).*x(2)).*(L2.*m2.*cos(x(1)+(-1).*x(2)).*((-1).*u(1)+ ...
  u(2)+L1.*(g.*m1.*sin(x(1))+m2.*(g.*sin(x(1))+L2.*x(4).^2.*sin(x(1)+( ...
  -1).*x(2)))))+L1.*(m1+m2).*(u(2)+(-1).*L2.*m2.*((-1).*L1.*x(3).^2.* ...
  sin(x(1)+(-1).*x(2))+g.*sin(x(2))))))+m2.*cos((1/2).*dt.*(x(3)+( ...
  -1).*x(4))+x(1)+(-1).*x(2)).*(16.*L1.^3.*L2.^2.*(m1+m2.*sin(x(1)+( ...
  -1).*x(2)).^2).^2+(2.*dt.*L2.^2.*m2+2.*dt.*L1.^2.*(m1+m2)+4.*dt.* ...
  L1.*L2.*m2.*cos(x(1)+(-1).*x(2))).*(2.*dt.*L2.*((-1).*u(1)+u(2))+2.* ...
  dt.*L1.*(u(2).*cos(x(1)+(-1).*x(2))+L2.^2.*m2.*x(4).^2.*sin(x(1)+( ...
  -1).*x(2))+g.*L2.*(m1.*sin(x(1))+m2.*cos(x(2)).*sin(x(1)+(-1).* ...
  x(2))))+L1.^2.*L2.*x(3).*((-4).*m1+m2.*((-4).*sin(x(1)+(-1).*x(2)) ...
  .^2+dt.*x(3).*sin(2.*(x(1)+(-1).*x(2)))))).*sin((1/2).*dt.*(x(3)+( ...
  -1).*x(4))+x(1)+(-1).*x(2))+(2.*dt.*L2+2.*dt.*L1.*cos(x(1)+(-1).* ...
  x(2))).*(2.*dt.*L2.^2.*m2.*((-1).*u(1)+u(2))+2.*dt.*L1.^3.*L2.*m2.*( ...
  m1+m2).*x(3).^2.*sin(x(1)+(-1).*x(2))+dt.*L1.*L2.*m2.*((-2).*(u(1)+( ...
  -2).*u(2)).*cos(x(1)+(-1).*x(2))+2.*L2.^2.*m2.*x(4).^2.*sin(x(1)+( ...
  -1).*x(2))+2.*g.*L2.*(m1.*sin(x(1))+m2.*cos(x(2)).*sin(x(1)+(-1).* ...
  x(2))))+L1.^2.*(2.*dt.*(m1+m2).*u(2)+2.*dt.*g.*L2.*m2.*(m1+m2).*cos( ...
  x(1)).*sin(x(1)+(-1).*x(2))+L2.^2.*m2.*((-4).*m1.*(x(3)+(-1).* ...
  x(4))+m2.*((-4).*x(3).*sin(x(1)+(-1).*x(2)).^2+4.*x(4).*sin(x(1)+( ...
  -1).*x(2)).^2+dt.*x(3).^2.*sin(2.*(x(1)+(-1).*x(2)))+dt.*x(4).^2.* ...
  sin(2.*(x(1)+(-1).*x(2))))))).*sin((1/2).*dt.*(x(3)+(-1).*x(4))+ ...
  x(1)+(-1).*x(2))+(-2).*L1.*((dt.*L1.*(m1+m2)+dt.*L2.*m2.*cos(x(1)+ ...
  (-1).*x(2))).*(2.*dt.*L2.*((-1).*u(1)+u(2))+2.*dt.*L1.*(u(2).*cos(x(1)+( ...
  -1).*x(2))+L2.^2.*m2.*x(4).^2.*sin(x(1)+(-1).*x(2))+g.*L2.*(m1.* ...
  sin(x(1))+m2.*cos(x(2)).*sin(x(1)+(-1).*x(2))))+L1.^2.*L2.*x(3).*( ...
  (-4).*m1+m2.*((-4).*sin(x(1)+(-1).*x(2)).^2+dt.*x(3).*sin(2.*( ...
  x(1)+(-1).*x(2)))))).*sin((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+(-1).* ...
  x(2))+(2.*dt.*L2+2.*dt.*L1.*cos(x(1)+(-1).*x(2))).*sin((1/2).*dt.* ...
  (x(3)+(-1).*x(4))+x(1)+(-1).*x(2)).*(2.*L1.*L2.^2.*m2.*x(4).*(m1+ ...
  m2.*sin(x(1)+(-1).*x(2)).^2)+dt.*L2.*m2.*cos(x(1)+(-1).*x(2)).*(( ...
  -1).*u(1)+u(2)+L1.*(g.*m1.*sin(x(1))+m2.*(g.*sin(x(1))+L2.*x(4).^2.* ...
  sin(x(1)+(-1).*x(2)))))+dt.*L1.*(m1+m2).*(u(2)+(-1).*L2.*m2.*((-1).* ...
  L1.*x(3).^2.*sin(x(1)+(-1).*x(2))+g.*sin(x(2))))))));

b41=(1/16).*dt.*L1.^(-3).*L2.^(-4).*m2.^(-1).*(m1+m2.*sin(x(1)+(-1).* ...
  x(2)).^2).^(-2).*(m1+m2.*sin((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+( ...
  -1).*x(2)).^2).^(-1).*(4.*L1.*cos((1/2).*dt.*(x(3)+(-1).*x(4))+ ...
  x(1)+(-1).*x(2)).*((-4).*L1.*L2.^3.*m2.*(m1+m2.*sin(x(1)+(-1).* ...
  x(2)).^2).^2+(-4).*dt.*L1.*L2.^3.*m2.^2.*x(4).*cos(x(1)+(-1).* ...
  x(2)).*(m1+m2.*sin(x(1)+(-1).*x(2)).^2).*sin((1/2).*dt.*(x(3)+(-1) ...
  .*x(4))+x(1)+(-1).*x(2))+(-2).*dt.^2.*L2.*m2.*cos(x(1)+(-1).*x(2)) ...
  .*sin((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+(-1).*x(2)).*(L2.*m2.*cos( ...
  x(1)+(-1).*x(2)).*((-1).*u(1)+u(2)+L1.*(g.*m1.*sin(x(1))+m2.*(g.*sin( ...
  x(1))+L2.*x(4).^2.*sin(x(1)+(-1).*x(2)))))+L1.*(m1+m2).*(u(2)+(-1).* ...
  L2.*m2.*((-1).*L1.*x(3).^2.*sin(x(1)+(-1).*x(2))+g.*sin(x(2))))))+ ...
  (m1+m2).*(((-2).*dt.*L2.^2.*m2+(-2).*dt.*L1.*L2.*m2.*cos(x(1)+(-1) ...
  .*x(2))).*(2.*dt.*L2.*((-1).*u(1)+u(2))+2.*dt.*L1.*(u(2).*cos(x(1)+(-1) ...
  .*x(2))+L2.^2.*m2.*x(4).^2.*sin(x(1)+(-1).*x(2))+g.*L2.*(m1.*sin( ...
  x(1))+m2.*cos(x(2)).*sin(x(1)+(-1).*x(2))))+L1.^2.*L2.*x(3).*((-4) ...
  .*m1+m2.*((-4).*sin(x(1)+(-1).*x(2)).^2+dt.*x(3).*sin(2.*(x(1)+( ...
  -1).*x(2)))))).*sin((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+(-1).*x(2))+ ...
  (-2).*dt.*L2.*(2.*dt.*L2.^2.*m2.*((-1).*u(1)+u(2))+2.*dt.*L1.^3.*L2.* ...
  m2.*(m1+m2).*x(3).^2.*sin(x(1)+(-1).*x(2))+dt.*L1.*L2.*m2.*((-2).* ...
  (u(1)+(-2).*u(2)).*cos(x(1)+(-1).*x(2))+2.*L2.^2.*m2.*x(4).^2.*sin( ...
  x(1)+(-1).*x(2))+2.*g.*L2.*(m1.*sin(x(1))+m2.*cos(x(2)).*sin(x(1)+ ...
  (-1).*x(2))))+L1.^2.*(2.*dt.*(m1+m2).*u(2)+2.*dt.*g.*L2.*m2.*(m1+m2) ...
  .*cos(x(1)).*sin(x(1)+(-1).*x(2))+L2.^2.*m2.*((-4).*m1.*(x(3)+(-1) ...
  .*x(4))+m2.*((-4).*x(3).*sin(x(1)+(-1).*x(2)).^2+4.*x(4).*sin( ...
  x(1)+(-1).*x(2)).^2+dt.*x(3).^2.*sin(2.*(x(1)+(-1).*x(2)))+dt.* ...
  x(4).^2.*sin(2.*(x(1)+(-1).*x(2))))))).*sin((1/2).*dt.*(x(3)+(-1) ...
  .*x(4))+x(1)+(-1).*x(2))+(-2).*L1.*((-1).*dt.*L2.*m2.*cos(x(1)+( ...
  -1).*x(2)).*(2.*dt.*L2.*((-1).*u(1)+u(2))+2.*dt.*L1.*(u(2).*cos(x(1)+( ...
  -1).*x(2))+L2.^2.*m2.*x(4).^2.*sin(x(1)+(-1).*x(2))+g.*L2.*(m1.* ...
  sin(x(1))+m2.*cos(x(2)).*sin(x(1)+(-1).*x(2))))+L1.^2.*L2.*x(3).*( ...
  (-4).*m1+m2.*((-4).*sin(x(1)+(-1).*x(2)).^2+dt.*x(3).*sin(2.*( ...
  x(1)+(-1).*x(2)))))).*sin((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+(-1).* ...
  x(2))+(-2).*dt.*L2.*sin((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+(-1).* ...
  x(2)).*(2.*L1.*L2.^2.*m2.*x(4).*(m1+m2.*sin(x(1)+(-1).*x(2)).^2)+ ...
  dt.*L2.*m2.*cos(x(1)+(-1).*x(2)).*((-1).*u(1)+u(2)+L1.*(g.*m1.*sin( ...
  x(1))+m2.*(g.*sin(x(1))+L2.*x(4).^2.*sin(x(1)+(-1).*x(2)))))+dt.* ...
  L1.*(m1+m2).*(u(2)+(-1).*L2.*m2.*((-1).*L1.*x(3).^2.*sin(x(1)+(-1).* ...
  x(2))+g.*sin(x(2))))))));


b42=(1/16).*dt.*L1.^(-3).*L2.^(-4).*m2.^(-1).*(m1+m2.*sin(x(1)+(-1).* ...
  x(2)).^2).^(-2).*(m1+m2.*sin((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+( ...
  -1).*x(2)).^2).^(-1).*(4.*L1.*cos((1/2).*dt.*(x(3)+(-1).*x(4))+ ...
  x(1)+(-1).*x(2)).*(4.*L1.*L2.^3.*m2.*(m1+m2.*sin(x(1)+(-1).*x(2)) ...
  .^2).^2+4.*dt.*L1.*L2.^2.*m2.*x(4).*(L1.*(m1+m2)+L2.*m2.*cos(x(1)+ ...
  (-1).*x(2))).*(m1+m2.*sin(x(1)+(-1).*x(2)).^2).*sin((1/2).*dt.*( ...
  x(3)+(-1).*x(4))+x(1)+(-1).*x(2))+2.*dt.^2.*(L1.*(m1+m2)+L2.*m2.* ...
  cos(x(1)+(-1).*x(2))).*sin((1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+(-1) ...
  .*x(2)).*(L2.*m2.*cos(x(1)+(-1).*x(2)).*((-1).*u(1)+u(2)+L1.*(g.*m1.* ...
  sin(x(1))+m2.*(g.*sin(x(1))+L2.*x(4).^2.*sin(x(1)+(-1).*x(2)))))+ ...
  L1.*(m1+m2).*(u(2)+(-1).*L2.*m2.*((-1).*L1.*x(3).^2.*sin(x(1)+(-1).* ...
  x(2))+g.*sin(x(2))))))+(m1+m2).*(16.*L1.^3.*L2.^2.*(m1+m2.*sin( ...
  x(1)+(-1).*x(2)).^2).^2+(2.*dt.*L2.^2.*m2+2.*dt.*L1.^2.*(m1+m2)+ ...
  4.*dt.*L1.*L2.*m2.*cos(x(1)+(-1).*x(2))).*(2.*dt.*L2.*((-1).*u(1)+ ...
  u(2))+2.*dt.*L1.*(u(2).*cos(x(1)+(-1).*x(2))+L2.^2.*m2.*x(4).^2.*sin( ...
  x(1)+(-1).*x(2))+g.*L2.*(m1.*sin(x(1))+m2.*cos(x(2)).*sin(x(1)+( ...
  -1).*x(2))))+L1.^2.*L2.*x(3).*((-4).*m1+m2.*((-4).*sin(x(1)+(-1).* ...
  x(2)).^2+dt.*x(3).*sin(2.*(x(1)+(-1).*x(2)))))).*sin((1/2).*dt.*( ...
  x(3)+(-1).*x(4))+x(1)+(-1).*x(2))+(2.*dt.*L2+2.*dt.*L1.*cos(x(1)+( ...
  -1).*x(2))).*(2.*dt.*L2.^2.*m2.*((-1).*u(1)+u(2))+2.*dt.*L1.^3.*L2.* ...
  m2.*(m1+m2).*x(3).^2.*sin(x(1)+(-1).*x(2))+dt.*L1.*L2.*m2.*((-2).* ...
  (u(1)+(-2).*u(2)).*cos(x(1)+(-1).*x(2))+2.*L2.^2.*m2.*x(4).^2.*sin( ...
  x(1)+(-1).*x(2))+2.*g.*L2.*(m1.*sin(x(1))+m2.*cos(x(2)).*sin(x(1)+ ...
  (-1).*x(2))))+L1.^2.*(2.*dt.*(m1+m2).*u(2)+2.*dt.*g.*L2.*m2.*(m1+m2) ...
  .*cos(x(1)).*sin(x(1)+(-1).*x(2))+L2.^2.*m2.*((-4).*m1.*(x(3)+(-1) ...
  .*x(4))+m2.*((-4).*x(3).*sin(x(1)+(-1).*x(2)).^2+4.*x(4).*sin( ...
  x(1)+(-1).*x(2)).^2+dt.*x(3).^2.*sin(2.*(x(1)+(-1).*x(2)))+dt.* ...
  x(4).^2.*sin(2.*(x(1)+(-1).*x(2))))))).*sin((1/2).*dt.*(x(3)+(-1) ...
  .*x(4))+x(1)+(-1).*x(2))+(-2).*L1.*((dt.*L1.*(m1+m2)+dt.*L2.*m2.* ...
  cos(x(1)+(-1).*x(2))).*(2.*dt.*L2.*((-1).*u(1)+u(2))+2.*dt.*L1.*(u(2).* ...
  cos(x(1)+(-1).*x(2))+L2.^2.*m2.*x(4).^2.*sin(x(1)+(-1).*x(2))+g.* ...
  L2.*(m1.*sin(x(1))+m2.*cos(x(2)).*sin(x(1)+(-1).*x(2))))+L1.^2.* ...
  L2.*x(3).*((-4).*m1+m2.*((-4).*sin(x(1)+(-1).*x(2)).^2+dt.*x(3).* ...
  sin(2.*(x(1)+(-1).*x(2)))))).*sin((1/2).*dt.*(x(3)+(-1).*x(4))+ ...
  x(1)+(-1).*x(2))+(2.*dt.*L2+2.*dt.*L1.*cos(x(1)+(-1).*x(2))).*sin( ...
  (1/2).*dt.*(x(3)+(-1).*x(4))+x(1)+(-1).*x(2)).*(2.*L1.*L2.^2.*m2.* ...
  x(4).*(m1+m2.*sin(x(1)+(-1).*x(2)).^2)+dt.*L2.*m2.*cos(x(1)+(-1).* ...
  x(2)).*((-1).*u(1)+u(2)+L1.*(g.*m1.*sin(x(1))+m2.*(g.*sin(x(1))+L2.* ...
  x(4).^2.*sin(x(1)+(-1).*x(2)))))+dt.*L1.*(m1+m2).*(u(2)+(-1).*L2.* ...
  m2.*((-1).*L1.*x(3).^2.*sin(x(1)+(-1).*x(2))+g.*sin(x(2))))))));

B= [b11,b12;b21,b22;b31,b32;b41,b42];