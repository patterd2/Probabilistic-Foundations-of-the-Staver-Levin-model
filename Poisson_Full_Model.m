% Matrix simulation of Poisson Forest-Grass Markov Chain
clear all;
tic

%rng(12352357);

%% parameter values
mu = 0.1;
nu = 0.05;
alpha = 0.3;
beta = 0.4;

w0 = 0.9;
w1 = 0.4;
s1 = 0.01;
t1 = 0.4;

t2 = 0.4;
f0 = 0.1;
f1 = 0.9;
s2 = 0.05;
%% setup, function definitions, etc.
sites = 1500;
P=1;                % Number of Patches
N=sites*ones(1,P);    % Number of Sites / Patch
NTot=sum(N);        % Total number of sites
T=300;
dt=0.01;
t0=0;
W=1;
J_F = alpha;
J_T = beta;
phi= @(x) (f0+(f1-f0)./(1+exp(-(x-t2)/s2)));
omega = @(x) (w0 + (w1-w0)./(1 + exp(-(x-t1)/s1)));
Ntimes = length(t0:dt:(T-dt));
% need to decide the relative proportions of cover to start with:
p0_grass = 0.25;
p0_sap = 0.25;
p0_savanna = 0.25;
%% Gillespie algorithm
Times=[0];
Grass = zeros(NTot,1);
Saplings = zeros(NTot,1);
Savanna = zeros(NTot,1);
Forest = zeros(NTot,1);

rand_initial = rand(N(1),1);
% allocate the initial cover according to prescribed proportions
Grass(:,1) = rand_initial < p0_grass;
Saplings(:,1) = (rand_initial < p0_grass + p0_sap) & (rand_initial > p0_grass);
Savanna(:,1) = (rand_initial < p0_grass + p0_sap + p0_savanna) & (rand_initial > p0_sap + p0_grass);
Forest(:,1) = rand_initial > p0_grass + p0_sap + p0_savanna;

k=0;
t=0;
% work out the transition intensities, want to keep track of these to
% plot at the end as well for mean-field comparison
grass_prop(1,k+1) = mean(Grass(:,1));
sap_prop(1,k+1) = mean(Saplings(:,1));
savanna_prop(1,k+1) = mean(Savanna(:,1));
forest_prop(1,k+1) = mean(Forest(:,1));
while (t<T)
    k=k+1;
    % compute the total transition intensity across the domain
    totalIntensity=sum(J_F*(forest_prop(1,k))*(1-Forest(:,k)) ... % not forest to forest
        + phi(W*grass_prop(1,k))*Forest(:,k) ... % forest to grass
        + J_T*savanna_prop(1,k).*Grass(:,k) ... % grass to sapling
        + omega(W*grass_prop(1,k))*Saplings(:,k) ... % sapling to savanna tree
        + mu*Saplings(:,k) + nu*Savanna(:,k)); % saplings/savanna to grass via natural mortality
    NextEvent=-log(1-rand())/totalIntensity;
    t=t+NextEvent; % compute the next time at which a site will change state
    Times(end+1)=t;
    if t>T
        Times(end) = T;
        grass_prop(1,k+1) = mean(Grass(:,k));
        sap_prop(1,k+1) = mean(Saplings(:,k));
        savanna_prop(1,k+1) = mean(Savanna(:,k));
        forest_prop(1,k+1) = mean(Forest(:,k));
        break;
    end
    % Is it a forest site which changes?
    norm_factor_switch = forest_prop(1,k)*phi(W*grass_prop(1,k))...
        + (1-forest_prop(1,k))*J_F*forest_prop(1,k) + J_T*savanna_prop(1,k)*grass_prop(1,k)...
        + omega(W*grass_prop(1,k))*sap_prop(1,k) + mu*sap_prop(1,k) + nu*savanna_prop(1,k);
    temp_rand = rand()*norm_factor_switch; % rather than normalise each time, just multiply once
    % 4 switches are now possible, need to select one with appropriate
    % probabilities
    ForestSwitch = temp_rand < (forest_prop(1,k)*phi(W*grass_prop(1,k)));
    SavannaSwitch = temp_rand < (forest_prop(1,k)*phi(W*grass_prop(1,k))+nu*savanna_prop(1,k)+J_F*savanna_prop(1,k)*forest_prop(1,k))...
        && temp_rand > (forest_prop(1,k)*phi(W*grass_prop(1,k)));
    GrassSwitch = temp_rand < (forest_prop(1,k)*phi(W*grass_prop(1,k))+nu*savanna_prop(1,k)+J_F*savanna_prop(1,k)*forest_prop(1,k)+grass_prop(1,k)*J_F*forest_prop(1,k)...
        +J_T*savanna_prop(1,k)*grass_prop(1,k))...
        && temp_rand > forest_prop(1,k)*phi(W*grass_prop(1,k))+nu*savanna_prop(1,k)+J_F*savanna_prop(1,k)*forest_prop(1,k);
    if ForestSwitch
        FF = find(Forest(:,end)==1); % find the forest sites
        Forest(:,k+1) = Forest(:,k);
        temp_switch = FF(randi(length(FF),1)); % pick one to switch
        Forest(temp_switch,k+1) = 0; % switch the selected forest site to grass
        Grass(:,k+1) = Grass(:,k);
        Grass(temp_switch,k+1) = 1; % update the grass vector to reflect the switch
        Saplings(:,k+1) = Saplings(:,k); % saplings and savanna's are unchanged
        Savanna(:,k+1) = Savanna(:,k);
    end
    if SavannaSwitch
        TT = find(Savanna(:,end)==1); % find the savanna sites
        Savanna(:,k+1) = Savanna(:,k);
        temp_switch = TT(randi(length(TT),1)); % pick one to switch
        Savanna(temp_switch,k+1) = 0;
        % Will it go to grass (mortality) or forest (seed invasion)?
        if rand() < nu*savanna_prop(1,k)/(nu*savanna_prop(1,k) + J_F*savanna_prop(1,k)*forest_prop(1,k)) % then it's a switch to grass
            Grass(:,k+1) = Grass(:,k);
            Saplings(:,k+1) = Saplings(:,k);
            Forest(:,k+1) = Forest(:,k);
            Grass(temp_switch,k+1) = 1;
        else % if not, it's a switch to forest
            Grass(:,k+1) = Grass(:,k);
            Saplings(:,k+1) = Saplings(:,k);
            Forest(:,k+1) = Forest(:,k);
            Forest(temp_switch,k+1) = 1;
        end
    end
    if GrassSwitch
        GG = find(Grass(:,end)==1); % find the grass sites
        Grass(:,k+1) = Grass(:,k);
        temp_switch = GG(randi(length(GG),1)); % pick one to switch
        Grass(temp_switch,k+1) = 0;
        % Will it go to saplings or forest?
        if rand() < grass_prop(1,k)*J_F*forest_prop(1,k)/(grass_prop(1,k)*J_F*forest_prop(1,k)+J_T*savanna_prop(1,k)*grass_prop(1,k)) % then it's a switch to forest
            Savanna(:,k+1) = Savanna(:,k);
            Saplings(:,k+1) = Saplings(:,k);
            Forest(:,k+1) = Forest(:,k);
            Forest(temp_switch,k+1) = 1;
        else % if not, it's a switch to saplings
            Savanna(:,k+1) = Savanna(:,k);
            Saplings(:,k+1) = Saplings(:,k);
            Forest(:,k+1) = Forest(:,k);
            Saplings(temp_switch,k+1) = 1;
        end
    end
    if ForestSwitch+SavannaSwitch+GrassSwitch==0 && sum(Saplings(:,end)) > 0  % then a sapling is switching to either savanna, forest or grass
        SS = find(Saplings(:,end)==1); % find the sapling sites
        Saplings(:,k+1) = Saplings(:,k);
        temp_switch = SS(randi(length(SS),1)); % pick one to switch
        Saplings(temp_switch,k+1) = 0;
        % Will it go to grass, forest or a savanna tree?
        temp_rand2 = rand();
        if temp_rand2 < mu*sap_prop(1,k)/(mu*sap_prop(1,k)+omega(W*grass_prop(1,k))*sap_prop(1,k)+J_F*forest_prop(1,k)*sap_prop(1,k)) % then it's a switch to grass by natural mortality
            Grass(:,k+1) = Grass(:,k);
            Savanna(:,k+1) = Savanna(:,k);
            Forest(:,k+1) = Forest(:,k);
            Grass(temp_switch,k+1) = 1;
        elseif temp_rand2 < (mu*sap_prop(1,k)+J_F*forest_prop(1,k)*sap_prop(1,k))/(mu*sap_prop(1,k)+omega(W*grass_prop(1,k))*sap_prop(1,k)+J_F*forest_prop(1,k)*sap_prop(1,k))...
                && temp_rand2 > mu*sap_prop(1,k)/(mu*sap_prop(1,k)+omega(W*grass_prop(1,k))*sap_prop(1,k)+J_F*forest_prop(1,k)*sap_prop(1,k)) % then it goes to forest
            Grass(:,k+1) = Grass(:,k);
            Savanna(:,k+1) = Savanna(:,k);
            Forest(:,k+1) = Forest(:,k);
            Forest(temp_switch,k+1) = 1;
        else % if not, it's a switch to savanna
            Forest(:,k+1) = Forest(:,k);
            Grass(:,k+1) = Grass(:,k);
            Savanna(:,k+1) = Savanna(:,k);
            Savanna(temp_switch,k+1) = 1;
        end
    end
    % update the cover proportions for next loop
    grass_prop(1,k+1) = mean(Grass(:,k+1));
    sap_prop(1,k+1) = mean(Saplings(:,k+1));
    savanna_prop(1,k+1) = mean(Savanna(:,k+1));
    forest_prop(1,k+1) = mean(Forest(:,k+1));
end
toc

grass_new = interp1(Times,grass_prop,0:dt:Times(end));
sap_new = interp1(Times,sap_prop,0:dt:Times(end));
savanna_new = interp1(Times,savanna_prop,0:dt:Times(end));
forest_new = interp1(Times,forest_prop,0:dt:Times(end));
%%
figure(2);
plot(0:dt:Times(end), grass_new);
hold on;
plot(0:dt:Times(end), sap_new);
hold on;
plot(0:dt:Times(end), savanna_new);
hold on;
plot(0:dt:Times(end), forest_new);
ylim([0 1]);
legend('Grass','Saplings','Savanna','Forest');

%%
figure(3);
imagesc(Grass(1:1000,:)+2*Forest(1:1000,:)+3*Savanna(1:1000,:)+4*Saplings(1:1000,:))

