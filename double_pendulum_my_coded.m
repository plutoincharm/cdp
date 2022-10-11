clear; clc; close all

%n = 8; % dimensions of system
global m1 m2 L1 L2 g val Nx Nu pert
m1=50.172;
m2=7.4;
m2=3.411;
L2=0.4418;L3=0.4033;
r2=0.1872;r3=0.1738;
g = 9.81; % gravity
pert = 0.0001;

Nx = 8;
Nu  = 4;
Tf = 1;
dt = 0.01;
Nt = round(Tf/dt)+1;



dynamics_midpoint = @(x,u,dt) x + fun_xdot(x + fun_xdot(x,u)*dt/2,u)*dt;
%{
A_midpoint = @(x,dt) [(1 + g*cos(x(1))*(dt^2)/(2*l)) (dt - mu*(dt^2)/(2*m*l^2));
                      (g*cos(x(1) + x(2)*dt/2)*dt/l - mu*g*cos(x(1))*(dt^2)/(2*m*l^3)) (1 + g*cos(x(1) + x(2)*dt/2)*(dt^2)/(2*l) - mu*dt/(m*l^2) + (mu^2)*(dt^2)/(2*(m^2)*l^4))];

B_midpoint = @(x,dt) [(dt^2)/(2*m*l^2); 
                      (-mu*(dt^2)/(2*(m^2)*l^4) + dt/(m*l^2))];

%}




% initial conditions
x1 = 0.5; y1 = 0.5; tht2=deg2rad(20);tht3=deg2rad(35);vx1=1.5;vy1=1.5;omg2=0.5;omg3=0.5;
x0 = [x1;y1;tht2;tht3;vx1;vy1;omg2;omg3];

% goal
xf = [2;2;pi;pi/2;2;2;3;3]; 

% costs
%Q = 1e-5*eye(4);
Q =  1e-5*eye(Nx);
Qf = 15*eye(Nx);
R = 5*1e-4*eye(Nu);
I = eye(4);

e_dJ = 1e-12;

% initialization
u = zeros(Nu,Nt-1);
x = zeros(Nx,Nt);
x_prev = zeros(Nx,Nt);
x(:,1) = x0;

% A = fun_amat(x(:,5),u(:,5),dt)
% B = fun_bmat(x(:,5),u(:,5),dt)
%  pause()
% first roll-out
for k = 2:Nt-1
        x(:,k) = dynamics_midpoint(x(:,k-1),u(:,k-1),dt);
        pause()
        
end
val = x;
tdt = 0.2:0.01:0.74;
tt = 0.23:dt:1.58;
fr = 1:1:125;
cr = 1:1:101;
hr = 1:1:101;


%%% theta2
figure;
plot(hr,val(3,1:101),'b-','LineWidth',1);
grid on;

xlabel('time (sec) \rightarrow');
ylabel('theta2 (N) \rightarrow');
legend('calc','dataset');


%%% theta3
figure;
plot(hr,val(4,1:101),'b-','LineWidth',1);
grid on;

xlabel('time (sec) \rightarrow');
ylabel('theta3 (N) \rightarrow');
legend('calc','dataset');















% original cost
J = 0;
for k = 1:Nt-1
    J = J + (x(:,k)-xf)'*Q*( x(:,k)-xf) + u(:,k)'*R*u(:,k);
end
disp('Original cost:')
J = 0.5*(J + (x(:,Nt)-xf)'*Qf*(x(:,Nt)-xf))

%{

pause()

%%%%%%%%%%%%%%%% ILQR Algorithm  %%%%%%%%%%%%%%%%%%%%
p = ones(Nx,Nt);
P = zeros(Nx,Nx,Nt);
%d = ones(Nu,Nu,Nt-1);
d = ones(Nu,Nt-1);
K = zeros(Nu,Nx,Nt-1);
%pdim = ones(Nx,Nu,Nt);
dJ = 1.0;  % change in cost

xn = zeros(Nx,Nt);
un = zeros(Nu,Nt-1);
% func g(dx,du) is perturbation of val func
% grad- g/ hessian-G of change in value fun
gx = zeros(Nx);
gu = zeros(Nu);
Gxx = zeros(Nx,Nx);
Guu = zeros(Nu,Nu);
Gxu = zeros(Nx,Nu);
Gux = zeros(Nu,Nx);

iter = 0;
while max(abs(d(:))) >  1e-3
    
    iter = iter +  1 

 %%%%% Backward Pass %%%%%
    dJ = 0.0;
    p(:,Nt) = Qf*(x(:,Nt)-xf);     %%% P is vx
    P(:,:,Nt) = Qf;                %%% P is vxx
    mu_reg = 0;
    for k = (Nt-1):-1:1
          %Calculate derivatives
           q = Q*( x(:,k)-xf);     % lx
           r = R*u(:,k);           % lu
            
            A = fun_amat(x(:,k), u(:,k),dt);
            B = fun_bmat(x(:,k), u(:,k),dt);

           %gradient of change in val fn
            gx = q + A'*p(:,k+1);% gx = dg/dx  
            gu = r + B'*p(:,k+1);% gu = dg/du
    
          %iLQR (Gauss-Newton) version
          %Hessian
             Gxx = Q + A'*(P(:,:,k+1))*A;
             Guu = R + B'*(P(:,:,k+1)+ mu_reg*eye(Nx))*B;
             Gxu = A'*P(:,:,k+1)*B;
             Gux = B'*(P(:,:,k+1) + mu_reg*eye(Nx))*A;     
             
             %beta = 0.1;
             log = issymmetric([Guu]);
             eigv = eig([Guu]);

          if any(eig(Guu)<0)
            mu_reg = mu_reg + 1;
            k = Nt-1;
            disp('regularized')
          end
        %{
              while (log==0) || all(eigv < 0) 
                    Gxx = Gxx + A'*beta*I*A
                    Guu = Guu + B'*beta*I*B
                    Gxu = Gxu + A'*beta*I*B
                    Gux = Gux + B'*beta*I*A
                    beta = 2*beta
                    %display("regularizing G")
                    display(beta)
                    log = issymmetric([Gxx Gxu; Gux Guu]);
                    eigv = eig([Gxx Gxu; Gux Guu]);
              end
         %}
            d(:,k) = Guu\gu;  % feedforward term
            K(:,:,k) = Guu\Gux; % feedback gain term
    
             p(:,k) = gx - K(:,:,k)'*gu + K(:,:,k)'*Guu*d(:,k) - Gxu*d(:,k);
             P(:,:,k) = Gxx + K(:,:,k)'*Guu*K(:,:,k) - Gxu*K(:,:,k) - K(:,:,k)'*Gux;
             dJ = dJ +  gu'*d(:,k);

          
    end
    disp('EOBP')
    pause()
  
%%%% End of Backward Pass %%%%%
      
    %Forward rollout with line search
    xn(:,1) = x(:,1);
    alpha = 1.0;
   for k = 1:(Nt-1)
        un(:,k) = u(:,k) - alpha*d(:,k) - (K(:,:,k)*(xn(:,k)-x(:,k)));
        xn(:,k+1) = dynamics_midpoint(xn(:,k),un(:,k),dt);
    end
    disp('EOFP')
    Jn = 0;
    for k = 1:Nt-1
        Jn = Jn + (xn(:,k)-xf)'*Q*(xn(:,k)-xf) + un(:,k)'*R*un(:,k);
    end
   Jn = 0.5*(Jn + (xn(:,Nt)-xf)'*Qf*(xn(:,Nt)-xf))
    
    
 
    while isnan(Jn) || Jn > (J - 1e-2*alpha*dJ)
        alpha = 0.5*alpha
        for k = 1:(Nt-1)
            un(:,k) = u(:,k) - alpha*d(:,k) - (K(:,:,k)*(xn(:,k)-x(:,k)));
            xn(:,k+1) = dynamics_midpoint(xn(:,k),un(:,k),dt);
           
        end
       
        Jn = 0;
        for k = 1:Nt-1
            Jn = Jn + (xn(:,k) - xf)'*Q*(xn(:,k) - xf) + un(:,k)'*R*un(:,k);
        end
     Jn = 0.5*(Jn + (xn(:,Nt) - xf)'*Qf*(xn(:,Nt) - xf))
    end
    
 
    J = Jn;
    x = xn;
    u = un;
  end
%}
   
    


















