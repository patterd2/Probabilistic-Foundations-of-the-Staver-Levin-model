% close all
tic
n=1000;
t=200;

M=zeros(n);

f0=0.1;
f1=0.9;
s2=0.05;
t2=0.4;

Phi = @(x) (f0+(f1-f0)./(1+exp(-(x-t2)/s2)));

% close all
Rmax=[];

alpha_vect=0:0.01:2;

Q_save=zeros(n-1,length(alpha_vect));
k=0;

for alpha=alpha_vect
    k=k+1;
    M=zeros(n);
    for i=2:(n-1)
        M(i,i+1)= alpha*(i)*(1-i/n);
        M(i,i-1)= (i)*Phi(1-i/n);
        M(i,i)=-(sum(M(i,:)));
    end
    M(n,n-1)=n*Phi(0);
    M(n,n)=-n*Phi(0);
%     for i=2:(n-1)
%         M(i,i+1)= alpha*i*(1-i/n);
%         M(i,i-1)= i*Phi(1-i/n);
%         M(i,i)=-(sum(M(i,:)));
%     end
%     M(n,n-1)=Phi(1);
%     M(n,n)=-Phi(1);
    Q=expm(t*M(2:end,2:end));
    %     Q=M(2:end,2:end);
    %     [V,D]=eigs((transpose(Q)),1);
    V=rand(size(Q,1),100);
    %V = ones(size(Q,1),1)/n;
    %     for i=1:ntimes
    V=transpose(Q)*V;
    V=mean(V,2);
    %     end
    Q_save(:,k)=V/sum(V);
end
%%
[AA,BB]=meshgrid(alpha_vect,1:(n-1));
figure(1);
surf(AA,BB,flipud(Q_save))
shading interp;
zlim([0 0.06])
caxis([0, 0.06])
colorbar;

figure(2);
h=pcolor(AA,BB,flipud(Q_save));
shading interp;
set(h, 'EdgeColor', 'none');
caxis([0, 0.02])
custom_map = [
    linspace(1,0,100)' linspace(1,0,100)' linspace(1,1,100)'];
colormap(custom_map);
colorbar;

toc
%figure(3);
%M1=sum(Q_save(1:ceil(0.4*n),:),1);
%M2=sum(Q_save(ceil(0.4*n)+1:end,:),1);
%plot(alpha_vect,M1);
%hold on
%plot(alpha_vect,M2);

%figure(3);
%plot(alpha_vect,sum(Q_save(1:floor(0.4*n),:),1));
%hold on
%plot(alpha_vect,sum(Q_save((floor(0.4*n)+1):end,:),1));
