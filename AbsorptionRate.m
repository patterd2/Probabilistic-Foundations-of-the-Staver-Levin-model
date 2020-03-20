%% The forest species goes extinct at an exponential rate, 
% plot the exponential extinction rate parameter as a function of the forest 
%birth rate parameter
tic;
n=1500;
t=1;
% close all
M=zeros(n);

f0=0.1;
f1=0.9;
s2=0.05;
t2=0.4;

Phi=@(x) (f0+(f1-f0)./(1+exp(-(x-t2)/s2)));

% close all
Rmax=[];
% close all
A = 0.1:0.01:0.49;
B = 0.5:0.0025:0.6;
C = 0.61:0.01:1;
alpha_vect=[A B C];
for alpha=alpha_vect
% alpha=1.5
    M=zeros(n); 
    for i=2:(n-1)
        M(i,i+1)= alpha*i*(1-i/n);
        M(i,i-1)= i*Phi(1-i/n);
        M(i,i)=-(sum(M(i,:)));
    end
    M(n,n-1)=Phi(0);
    M(n,n)=-Phi(0);

    E=eig(M(2:end,2:end));
    RE=E(abs(imag(E))<0.001);
    Rmax(end+1)=max(real(RE));
    
end
%% Plot
figure(1)
hold on
plot(alpha_vect,-Rmax);
hold on;
line([0.5466 0.5466], [-0.02 0.12]);

toc;
