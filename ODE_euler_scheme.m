%% Explicit Euler in time to solve the 4D ODE system
%% function definitions
omega = @(w_0, w_1, g, t_1, s_1) w_0 + (w_1-w_0)./(1 + exp(-(g-t_1)/s_1));
phi = @(f_0, f_1, g, t_2, s_2) f_0 + (f_1-f_0)./(1 + exp(-(g-t_2)/s_2));
%% model parameters (values from PNAS paper)
mu = 0.1;
nu = 0.05;
alpha = 0.25; 
beta = 0.4; 

w_0 = 0.9;
w_1 = 0.4;
t_1 = 0.4;
s_1 = 0.01;

f_0 = 0.1;
f_1 = 0.9;
t_2 = 0.4;
s_2 = 0.05;
%% numerical method parameters
h = 0.01; % time discretisation parameter
n = 50000; % number of time steps
tau = n*h; % simulations time domain is [0,tau]
%% initial conditions
G = 0.25*ones(1,n+1);
T = 0.25*ones(1,n+1);
S = 0.25*ones(1,n+1);
F = ones(1,n+1) - S - T - G;
%% The Euler scheme
alpha_0=0.1;
alpha_1=0.4;
Nalpha=10;
kalpha=0;
% Ntimes=length(0:dt:T);
close all

figure(1)
hold on
for alpha=linspace(alpha_0,alpha_1,Nalpha)
for i = 2:n+1

    G(1,i) = G(1,i-1) + h*( mu*S(1,i-1) + nu*T(1,i-1)...
        - alpha*F(1,i-1)*G(1,i-1) + ...
        phi(f_0, f_1, G(1,i-1), t_2, s_2).*F(1,i-1)...
        - beta*T(1,i-1)*G(1,i-1) );
    
    S(1,i) = S(1,i-1) + h*( -mu*S(1,i-1) + beta*T(1,i-1)*G(1,i-1)...
        - alpha*F(1,i-1)*S(1,i-1)...
        - omega(w_0, w_1, G(1,i-1), t_1, s_1).*S(1,i-1) );
    
    T(1,i) = T(1,i-1) + h*( - nu*T(1,i-1) ...
        + omega(w_0, w_1, G(1,i-1), t_1, s_1).*S(1,i-1)...
        -alpha*F(1,i-1)*T(1,i-1) );
    
    F(1,i) = 1 - S(1,i) - T(1,i) - G(1,i);
    
end
plot3(G,T,F);
% plot([0 1],[0 1]);
end
xlabel('Grass') 
ylabel('Savanna Tree') 
zlabel('Forest Tree') 

legend('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1')
%% plot the solution versus time
figure(2);
plot(0:h:tau,G);
hold on;
plot(0:h:tau,S);
hold on;
plot(0:h:tau,T);
hold on;
plot(0:h:tau,F);
ylim([0 1]);
legend('Grass','Saplings','Savanna','Forest');
