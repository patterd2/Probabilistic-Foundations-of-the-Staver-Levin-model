% Matrix simulation of Poisson Forest-Grass Markov Chain
clear all;
Times=[];
DO_PLOT=0;
t2=0.4;
f0=0.1;
f1=0.9;
s2=0.05;
% Current simulation is for a single patch, i.e. M = 1
sites=100;
alpha=0.05:0.05:2;
p0=0.05:0.025:0.975;
final_state = zeros(length(alpha),length(p0));
P=1;                % Number of Patches
N=sites*ones(1,P);    % Number of Sites / Patch
NTot=sum(N);        % Total number of sites
T=100;
dt=0.01;
t0=0;
MC=1;
W=1;
phi=@(x) (f0+(f1-f0)./(1+exp(-(x-t2)/s2)));
Ntimes=length(t0:dt:(T-dt));
tic
%% run enough simulations to compute the expected extinction time as a function of the number of sites
for ii = 1:length(alpha)
    J=alpha(ii);
    for kk = 1:length(p0)
        Times=[0];
        % This block of code runs one simulation via the Gillespie algorithm
        for iteration=1:MC
            %fprintf('Monte Carlo Simulation #=%d/%d\n',iteration,MC);
            Solution=zeros(NTot,1);
            Solution(:,1)=rand(N(1),1)>p0(kk);
            Fbar=mean(Solution(:,1));
            tStart=tic;
            k=0;
            t=0;
            while (t<T)
                k=k+1;
                totalIntensity=sum(J*Fbar*(1-Solution(:,k))+phi(W*(1-Fbar))*Solution(:,k));
                NextEvent=-log(1-rand())/totalIntensity;
                % decide the next time at which a site will change state
                t=t+NextEvent;
                Times(end+1)=t;
                % decide if it a forest or grass site which will flip
                ForestSwitch=rand()<Fbar*phi(W*(1-Fbar))/(Fbar*phi(W*(1-Fbar))+(1-Fbar)*J*Fbar);
                if ForestSwitch
                    FF=find(Solution(:,end)==1); % find all the forest sites
                    Solution(:,k+1)=Solution(:,k); % copy the solution vector
                    Solution(FF(randi(length(FF),1)),k+1)=0; % pick which forest site will flip
                else
                    % if it wasn't a forest site flipping, it was grass
                    GG=find(Solution(:,end)==0); % which sites are grass
                    Solution(:,k+1)=Solution(:,k);
                    Solution(GG(randi(length(GG),1)),k+1)=1; % pick a grass site to flip
                end
                Fbar=mean(Solution(:,k+1));
            end
        end
        final_state(ii,k) = mean(Solution(:,end));
    end
end
toc
for i = 1:length(alpha)
    hold on;
    scatter(alpha(i)*ones(1,length(final_state(1,:))),1-final_state(i,:),'*','b');
    hold on;
end
xlim([0 2]);
ylim([0 1]);
xlabel('\alpha','FontSize',20);
ylabel('Grass','FontSize',14);