

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% mom2pot_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test mom2pot.m

clear
% define two ellipsoids
nd=2;
fac=1;
mu1=fac*[1,5]'
L1=fac*[1 2;0 1]; 
C1=L1*L1'
mu2=fac*[-2,-4]'
L2=fac*[1 -2;0 1]; 
C2=L2*L2'

% create samples
N=10000;
data(:,:,1)=mu1*ones(1,N)+L1*randn(2,N);
data(:,:,2)=mu2*ones(1,N)+L2*randn(2,N);

% draw data and fac*sigma ellipse
figure(1);clf;hold on;
fac=2.5;
ellipse(mu1,fac,C1,'-b');
h=plot(data(1,:,1),data(2,:,1),'.b');
set(h,'markersize',3)
ellipse(mu2,fac,C2,'-r');
h=plot(data(1,:,2),data(2,:,2),'.r');
set(h,'markersize',3)


% create sufficient statistic
countx=[N,N];
sumx=zeros(nd,2);
Mx=zeros(nd,nd,2);
for t=1:2,
  sumx(:,t)=sum(data(:,:,t),2);
  Mx(:,:,t)=data(:,:,t)*data(:,:,t)';
end;

% create potential
reg=2;
[meanx,Rx,beta,Vmin]=mom2pot(countx,sumx,Mx,reg);


% estimate ellipse
C1e=inv(Rx(:,:,1)+Rx(:,:,1)')/beta(1)
mu1e=meanx(:,1)
C2e=inv(Rx(:,:,2)+Rx(:,:,2)')/beta(2)
mu2e=meanx(:,2)

% draw prediction
ellipse(mu1e,fac,C1e,'-.b');
ellipse(mu2e,fac,C2e,'--r');

% check shift
Vmin
0.5*log([det(2*pi*C1),det(2*pi*C2)])
0.5*log([det(2*pi*C1e),det(2*pi*C2e)])
