close all
L1=2;L2=3;
i = 1;
while i <  length(val(1,:))
  
tht1=val(1,i);tht2=val(1,i);
hx=0.5;hy=0.5;  



xk = hx  + L1*cos(tht1+pi); 
yk = hy  + L1*sin(tht1+pi); 
xa = hx + L1*cos(tht1+pi) + L2*cos(tht2+pi);
ya = hy + L1*sin(tht1+pi) + L2*sin(tht2+pi);


axis([-5 10 -5 10]) 
grid on
%axis([-0.5 1.5 -0.5 0.5]) 
%base =line([-1 2],[0.02 0.02],'LineWidth',1,'Color','black');
%pelvic =line([lhip_x(i) rhip_x(i)],[lhip_z(i) rhip_z(i)],'LineWidth',1,'Color','black');
  T=line([hx  xk],[hy  yk],'LineWidth',1,'Color','black');
  u=line([xk xa],[yk ya],'LineWidth',1,'Color','red');

  

 i = i +1;
 pause
end  



  