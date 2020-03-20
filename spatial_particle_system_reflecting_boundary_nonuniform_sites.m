%% Matrix simulation of Poisson Forest-Grass Markov Chain
clear all;

% RNG seed
%rng(1546793)

Times=[];
DO_PLOT=0;
No_Abs=1;

%% Spatial Parameters
L = 1; % working on [0,L]
a = 0.4; % minimal site density
b = 1.2; % site density slope parameter

%% dispersal and fire kernels
J_F_fun = @(x, a, sigma_F) exp( -( (a-x).^2)/(2*sigma_F^2) )/sqrt(2*pi*sigma_F^2);
W_fun = @(x, a, sigma_W) exp( -( (a-x).^2 )/(2*sigma_W^2) )/sqrt(2*pi*sigma_W^2);

sigma_F = 0.02; % seed dispersal radius forest trees
sigma_W = sigma_F; % fire radius

% Model Parameters
alpha=1.1;

t2=0.4;
f0=0.1;
f1=0.9;
s2=0.05;

sites=1000;

P=1;                % Number of Patches
N=sites*ones(1,P);    % Number of Sites / Patch
NTot=sum(N);        % Total number of sites

T=500;
dt=0.01;
t0=0;

MC=1;

J=alpha;
W=1;
phi=@(x) (f0+(f1-f0)./(1+exp(-(x-t2)/s2)));

Ntimes=length(t0:dt:(T-dt));

tic
nbins=100;
Histo=zeros(nbins,MC);
p_grass=0.5;
p_no_grass=0.5;
%% This block of code runs one simulation via the Gillespie algorithm
for iteration=1:MC
    Solution=zeros(NTot,1);
    Locations=sort( -(a/b)+sqrt( (a^2)/(b^2) + 2*rand(sites,1)/b)); % inverse transform method for random site locations
    % histogram(Locations, 50); % check the dsitribution of the sites
    [A,B]=meshgrid(Locations);
    [A_left,B_left]=meshgrid(Locations-L);
    [A_right,B_right]=meshgrid(Locations+L);
    J_Mat=(J_F_fun(A,B,sigma_F) + fliplr(J_F_fun(A,B_left,sigma_F)) + fliplr(J_F_fun(A,B_right,sigma_F)))/sites;
    W_Mat=(W_fun(A,B,sigma_W) + fliplr(W_fun(A,B_left,sigma_W)) + fliplr(W_fun(A,B_right,sigma_W)))/sites;
    
    % check that the normalization of the kernels is correct
    %[X,Y] = meshgrid(0:1/sites:L);
    %J_Mat_regular = J_F_fun(X,Y,sigma_F)/sites;
    %plot(max(sum(J_Mat_regular)))
    
    Solution(:,1)=rand(NTot,1)<p_grass;
    %No_grass=find((Locations>0.2*L).*(Locations<0.5*L));
    %N_no_grass=length(No_grass);
    %Solution(No_grass,1)=rand(N_no_grass,1)<p_no_grass;
    
    Times=[0];
    
    k=0;
    t=0;
    while (t<T)
        k=k+1;
        BirthRates=alpha*((J_Mat*(1-Solution(:,k)))).*Solution(:,k);
        DeathRates=phi(W_Mat*Solution(:,k)).*(1-Solution(:,k));
        %BirthRates=alpha*((J_Mat*(1-Solution(:,k)))).*Solution(:,k);
        %DeathRates=phi(W_Mat*Solution(:,k)).*(1-Solution(:,k));
        
        totalIntensity=sum(BirthRates+DeathRates);
        
        NextEvent=-log(1-rand())/totalIntensity;
        
        t=t+NextEvent;
        Times(end+1)=t;
        
        CDF=cumsum(BirthRates+DeathRates)/totalIntensity;
        U=rand();
        i=1;
        while U>CDF(i)
            i=i+1;
        end
        Solution(:,k+1)=Solution(:,k);
        Solution(i,k+1)=1-Solution(i,k);
    end
    k=k+1;
    Times(end)=T;
    %Solution(:,k+1)=Solution(:,k);
end
toc

%% Plots

% linearly interpolate solution values to put simulation on true timescale
sol_new = ones(sites, length(0:dt:Times(end)));
for i = 1:sites
    sol_new(i,:) = interp1(Times,Solution(i,:),0:dt:Times(end));
end

% need to put the spatial dimension on a linear grid as well
dx = 0.005; % spatial mesh
sol_new2 = ones(length(0:dx:Locations(end)),length(0:dt:Times(end)));
for i=1:length(0:dt:Times(end))
    sol_new2(:,i) = interp1(Locations,sol_new(:,i),0:dx:Locations(end));
end

% Space-time plot
figure(1);
imagesc(sol_new2);
custom_map = [1 1 1
    0 0.5 0];
colormap(custom_map);

% plot the average proportion of grass in the system 
% should be able to roughly observe the Maxwell point
figure(2);
plot(mean(sol_new2(:,floor(0.4*length(0:dt:Times(end))):end),2));
ylim([0 1]);
