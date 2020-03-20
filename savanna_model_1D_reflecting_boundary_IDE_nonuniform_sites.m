%% 1D Simulation of spatial savanna model with soil gradient (nonconvolution)
% Explicit Euler in time, 1D trapezoidal rule in space
% The boundaries are reflecting
tic
%% Numerical method parameters
L = 1; % working on [0,L]
N = 199; % N+1 grid points
delta = L/N;  % spatial discretization parameter
h = 0.05; % time discretisation parameter
n = 10001; % number of time steps
tau = (n-1)*h; % simulations time domain is [0,tau]
%% Function definitions
J_F_fun = @(x, y, sigma_F, a, b) (a + b*x).*exp( -( (y-x).^2)/(2*sigma_F^2))/sqrt(2*pi*sigma_F^2);
W_fun = @(x, y, sigma_W, a, b) (a + b*x).*exp( -( (y-x).^2 )/(2*sigma_W^2))/sqrt(2*pi*sigma_W^2) ;
phi = @(f_0_fun, f_1, g, t_2_fun, s_2) f_0_fun + (f_1-f_0_fun)./(1 + exp(-(g-t_2_fun)/s_2));
%% Model parameters
alpha = 1.1;

a = 0.4;
b = 1.2;

w_0_ref = 0.9;
w_1_ref = 0.4;
t_1_ref = 0.4;
s_1 = 0.01;
f_0_ref = 0.1;
f_1_ref = 0.9;
t_2_ref = 0.4;
s_2 = 0.05;% s_2 standard value is 0.05

sigma_F = 0.02; % seed dispersal radius forest trees
sigma_W = sigma_F; % fire spread radius
%% Set up the initial distributions of the cover types on the grid
% each row is one time step of the simulation
% solution is a block matrix [LB; SOL; RB] where LB and RB are fixed
G0 = 0.25+0.5*rand(1,N+1);
F0 = 0.25+0.5*rand(1,N+1);

CR = G0 + F0;
G0 = G0./CR;
F0 = F0./CR;
LB_G = G0(1,2:end);
RB_G = G0(1,1:end-1);
LB_F = F0(1,2:end);
RB_F = F0(1,1:end-1);

G = ones(n,3*N+1);
F = ones(n,3*N+1);

G(1,:) = [fliplr(LB_G) G0 fliplr(RB_G)];
F(1,:) = [fliplr(LB_F) F0 fliplr(RB_F)];

% compute the convolution for E
X = 0:delta:L;
X_L = X-L;
X_L = X_L(1,1:end-1);
X_R = X+L;
X_R = X_R(1,2:end);
E = ones(1,N+1);
temp_normalise = ones(1,N+1);
% Save tempW matrices to avoid computing them again
tempW = ones(N+1,3*N+1);
Trap = ones(1,3*N+1);
%Trap(1,3*N+1)=0.5;
%Trap(1,1)=0.5;
for i = 1:N+1
    tempW(i,:) =  delta*W_fun([X_L X X_R],(i-1)*delta,sigma_W,a,b);
    integrand = tempW(i,:).*(G(1,:));
    E(1,i) = sum(integrand.*Trap);
end
%% preallocaate some temp variables for efficiency
temp1 = ones(1,N+1);
tempF = ones(N+1,3*N+1);
%% pre-calculate the 4D convolution matrices
for k = 1:N+1
    tempF(k,:) = delta*J_F_fun([X_L X X_R],(k-1)*delta,sigma_F,a,b);
end
%% The numerical scheme
for i = 2:n
    % compute convolutions for this time step
    progressbar(i,n);
    for k = 1:N+1
        integrand1 = tempF(k,:).*F(i-1,:);
        temp1(1,k) = sum(integrand1.*Trap);
    end
    G(i,(N+1):2*N+1) = G(i-1,(N+1):2*N+1) + h*( ...
        - alpha.*temp1.*G(i-1,(N+1):2*N+1) + ...
        phi(f_0_ref, f_1_ref, E(i-1,:), t_2_ref, s_2).*F(i-1,(N+1):2*N+1));
    F(i,(N+1):2*N+1) = ones(1,N+1) - G(i,(N+1):2*N+1);
    % add in the boundary conditions
    G(i,1:N) = fliplr(G(i,N+2:2*N+1));
    G(i,2*N+2:3*N+1) = fliplr(G(i,(N+1):2*N));
    F(i,1:N) = fliplr(F(i,N+2:2*N+1));
    F(i,2*N+2:3*N+1) = fliplr(F(i,(N+1):2*N));
    for k = 1:N+1
        integrand = tempW(k,:).*G(i,:);
        E(i,k) = sum(integrand.*Trap);
    end
    % For numerical stability take max with zero and min with one
    G(i,:) = min(max(G(i,:),0),1);
    F(i,:) = min(max(F(i,:),0),1);
    E(i,:) = min(max(E(i,:),0),1);
end
% The following output is useful when trying to discern whether or not
% a solution is stationary in time
fprintf('\n');
fprintf('The maximum changes on the grid for each variable at the last time step were:\n');
fprintf(['G: ',num2str(max(abs(G(n,:)-G(n-1,:)))),'\n']);
fprintf(['F: ',num2str(max(abs(F(n,:)-F(n-1,:)))),'\n']);
toc
%% Space Time Plot
figure(3);
imagesc(G(:,N+1:2*N+1)');
shading interp;
custom_map = [
    linspace(1,0,100)' linspace(1,0.5,100)' linspace(1,0,100)'];
colormap(custom_map);
colorbar;

%% Cross section of the solution
figure(4)
plot(G(end,N+1:2*N+1));
ylim([0 1]);
xlim([0 length(G(end,N+1:2*N+1))]);
