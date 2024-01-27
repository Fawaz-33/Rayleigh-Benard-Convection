% Rayleigh Benard Convection
% Square Cavity
% Top surface -1 and bottom surface 1
% Simple algorithm (Finite volume method)

clear all
close all
clc

np = 51;                        %Increase grid points for greater Rayleigh Number
dom =1;
h = dom/(np-1);
x = 0:h:dom;
y = 0:h:dom;
Pr = 0.71;                      %Prandtl number
Ra = 1e4;                       %Rayleigh number
mu = sqrt(Pr/Ra);               
D = 1 / sqrt(Ra*Pr);        
T_top = -1;
T_bottom = 1;
% Under relaxation factors
alpha = 0.8;
alphaV = 0.9;
alphaU = 0.9;
alpha_p = 0.8;


%% Initial Condition

% Staggered grid
u(np+1,np) = 0;
u_star(np+1,np) = 0;
d_e(np+1,np) = 0;

v(np,np+1) = 0;
v_star(np,np+1) = 0;
d_n(np,np+1) = 0;

T(np,np) = 1;
T_old(np,np) = 1;

P(np+1,np+1) = 1;
P_star(np+1,np+1) = 1;
P_corr(np+1,np+1) = 1;

b(np+1,np+1) = 1;
u(1,:) = 0;

u_new(np+1,np) = 0;
v_new(np,np+1) = 0;
P_new(np+1,np+1) = 1;
u_new(1,:) = 0;


%% Governing Equations
error = 1;
iter = 0;
error_tol = 1e-6;
% Increase iteration for greater Ra number
while iter < 12000
    % X- Momentum
    for i = 2:np
      for  j = 2:(np-1)
          u_E = 0.5*(u(i,j)+u(i,j+1));
          u_W = 0.5*(u(i,j)+u(i,j-1));
          v_N = 0.5*(v(i-1,j)+v(i-1,j+1));
          v_S = 0.5*(v(i,j)+v(i,j+1));
          
          a_E = - (0.5*u_E*h) + mu;
          a_W =  (0.5*u_W*h) + mu;
          a_N = - (0.5*v_N*h) + mu;
          a_S =  (0.5*v_S*h) + mu;
          a_e = (0.5*u_E*h) - (0.5*u_W*h) + (0.5*v_N*h) - (0.5*v_S*h) + mu + mu + mu + mu;
          
          d_e(i,j) = (alphaU*h)/a_e;
          u_star(i,j) = alphaU/a_e * ( (a_E*u(i,j+1) + a_W*u(i,j-1) + a_N*u(i-1,j) + a_S*u(i+1,j)) + ((P(i,j+1)-P(i,j)) * h) ) + (1-alphaU)*u(i,j);
          
      end
    end
    
      % BC- X - Momentum
      u_star(1,:) =  - u_star(2,:); % Upper
      u_star(np+1,:) = -u_star(np,:); % Lower
      u_star(2:np,1) = 0; 
      u_star(2:np,np) = 0;
      
      % Y- Momentum
    for i = 2:(np-1)
      for  j = 2:np
          u_E = 0.5*(u(i,j)+u(i+1,j));
          u_W = 0.5*(u(i,j-1)+u(i+1,j-1));
          v_N = 0.5*(v(i,j)+v(i-1,j));
          v_S = 0.5*(v(i,j)+v(i+1,j));
          
          a_E = - (0.5*u_E*h) + mu;
          a_W =  (0.5*u_W*h) + mu;
          a_N = - (0.5*v_N*h) + mu;
          a_S =  (0.5*v_S*h) + mu;
          a_n = (0.5*u_E*h) - (0.5*u_W*h) + (0.5*v_N*h) - (0.5*v_S*h) + mu + mu + mu + mu;
          
          Boussinesq_sourceterm = 0.5*(T(i-1,j)+T(i,j))*h*h ;
          d_n(i,j) = h*alphaV/a_n; 
          v_star(i,j) = alphaV/a_n * ( (a_E*v(i,j+1) + a_W*v(i,j-1) + a_N*v(i-1,j) + a_S*v(i+1,j)) + ((P(i,j)-P(i+1,j))*h) + Boussinesq_sourceterm) + (1-alphaV)*v(i,j);
      end
    end
    
      % BC- Y - Momentum
      v_star(:,1) = -v_star(:,2); % Upper
      v_star(:,np+1) = -v_star(:,np); % Lower
      v_star(1,2:np) = 0;  
      v_star(np,2:np) = 0;
        
    
    %Zeroing pressure correction to start with
    P_corr(1:np+1,1:np+1) = 0;
    
    % Continuity equation
    for i = 2:np
        for j = 2:np
            a_E = - d_e(i,j)*h;
            a_W = - d_e(i,j-1)*h;
            a_N = - d_n(i-1,j)*h;
            a_S = - d_n(i,j)*h;
            a_P = a_E + a_W + a_N + a_S;
            
            b(i,j) = -(u_star(i,j)-u_star(i,j-1))*h + (v_star(i,j)-v_star(i-1,j))*h;
            P_corr(i,j) = ( a_E*P_corr(i,j+1) + a_W*P_corr(i,j-1) + a_N*P_corr(i-1,j) + a_S*P_corr(i+1,j) + b(i,j))/a_P;
            
        end
    end
    % Pressure field correction 
    for i = 2:np
        for j = 2:np
            P_new(i,j)= P(i,j) + alpha_p*P_corr(i,j);
                
        end
    end
    
    % BC - Continuity equation
    P_new(1,:)= P_new(2,:); % Upper
    P_new(np+1,:)= P_new(np,:); % Lower
    P_new(:,1)= P_new(:,2);
    P_new(:,np+1)= P_new(:,np);
    
    % Correction of velocities 
    
    for i = 2:np
      for  j = 2:(np-1)
          u_new(i,j) = u_star(i,j) + d_e(i,j)*(P_corr(i,j+1)-P_corr(i,j))*alpha;
      end
    end
    
    % X - BC
      u_new(1,:) =  - u_new(2,:); % Upper
      u_new(np+1,:) = -u_new(np,:); % Lower
      u_new(2:np,1) = 0; 
      u_new(2:np,np) = 0;
    
    for i = 2:np-1
      for  j = 2:np
          v_new(i,j) = v_star(i,j) + d_n(i,j)*(P_corr(i,j)-P_corr(i+1,j))*alpha;
      end
    end
    
    % Y - BC
      v_new(:,1) = -v_new(:,2); % Upper
      v_new(:,np+1) = -v_new(:,np); % Lower
      v_new(1,2:np) = 0;  
      v_new(np,2:np) = 0;
    
      % Temp
    for i = 2:np-1
      for  j = 2:np-1
          u_E = u_new(i,j+1);
          u_W = u_new(i,j);
          v_N = v_new(i+1,j);
          v_S = v_new(i,j);
          
          a_E = - (u_E*h) + D;
          a_W =  (u_W*h) + D;
          a_N = - (v_N*h) + D;
          a_S = (v_S*h) + D;
          a_T = -(u_E*h) + (u_W*h) - (v_N*h) + (v_S*h) + D + D + D + D;
          
          T(i,j) =  ( (a_E*T_old(i,j+1)+a_W*T_old(i,j-1)+a_N*T_old(i-1,j)+a_S*T_old(i+1,j)) )/a_T ;
      end
    end
    
    
      %%apply BCs
      T(1,:)= T_top;        %top wall
      T(np,:)=T_bottom;     %bottom wall
      T(:,1)=T(:,2);            
      T(:,np)=T(:,np-1);    
      
      
    % Error Calculation
    error = 0;
    for i = 2:(np - 1)
        for j =2:(np - 1)
           error = error + abs(b(i,j));
        end
    end
    
    % Assigning values
    u = u_new;
    v = v_new;
    P = P_new;
    T_old   = T;
    iter =iter+1;
       
end




%% Plotting

pcolor(x,y,T)
shading interp
colorbar
colormap('jet')
xlabel('x')
ylabel('y')
title('Rayleigh Benard Convection')
