%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% clusters_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test clusters.m
%


clear
clear mex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig=1;                  % figure number
new=1;                  % new problem?
cat=0;                  % number of categories per feature 
                        % (0 = numerical features)
d=2;                    % dimension of data
c=4;                    % number of classes
np=[2,1,2,1]*50;        % number of data points per class
var=[1 2 2 1]*10;       % std of classes (>10 needed for rounding!)
fac=[1 1 1 1]*200;    % approx box size 
                        % (fac/var=12 is interesting, >12 easier)

if 0,
  c=5;                    % number of classes
  np=[2,2,2,2,2]*50;      % number of data points per class
  var=[10 10 1 1 1]*10;   % std of classes (>10 needed for rounding!)
  fac=[2 2 1 1 1]*500;    % approx box size 
                          % (fac/var=12 is interesting, >12 easier)
end;

% set parameters for classifier
sup=0;                  % clustering is always without supervision
prt=1;                  % printlevel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cat==1,error('cat must be 0 or >1'); end;
draw=( d==2 & c<=35 & ~cat );
type ='o+x.*o+x.*o+x.*o+x.*o+x.*o+x.*o+x.*';
color='rrrrrbbbbbkkkkkmmmmmcccccyyyyyggggg';
color='rbkmcygrbkmcygrbkmcygrbkmcygrbkmcyg'; % preferable

save test_debug0 new sup draw

tic;
if new,

  % pick centers
  if d==1,
    centers=[1:c]';
  else
    centers=rand(c,d);
  end;

  % define numerical data
  % rounding data to integers is assumed to be be harmless!! 
  typ=zeros(1,d);
  nsample=0;
  for ic=1:c,
    ind=[nsample+1:nsample+np(ic)];
    data(ind,:)=fac(ic)*centers(ic+0*ind,:)+var(ic)*randn(np(ic),d);
    % targets(ind,ic)=-1+zeros(np(ic),1);      % unit costs
    targets(ind,1)=ic+zeros(np(ic),1);         % categorical data
    class(ind,1)=ic+0*ind';
    nsample=nsample+np(ic);
  end;

  if cat>0,
    % transform to categorical data
    typ=cat+typ;
    for i=1:d,
      data(:,i)=data(:,i)-min(data(:,i));
      data(:,i)=round(0.501+(cat-0.002)*(data(:,i)/max(data(:,i))));
    end;
  end;

  % pick (0,1)-weights
  weights=round(rand(1,size(targets,1)));

  save test_debug

else
  load test_debug
  load test_debug0

end;
disp([showtime,' for data generation']);

mag=1;
data=mag*data;


if draw,
  rad=2.5; % factor for ellipse drawing; best in [2,3]
  % plot data 
  figure(fig); clf; 
  subplot(1,2,1);hold on;
  for ic=1:c,
    ind=find(class==ic);
    plot(data(ind,1),data(ind,2),[type(ic),color(ic)]);
    if cat==0, 
      mu=fac(ic)*centers(ic,:);
      C=var(ic)^2*eye(d);
      ellipse(mu,rad,C,color(ic));
    end;
  end;
  xlim=[min(data(:,1)),max(data(:,1))];
  ylim=[min(data(:,2)),max(data(:,2))];
  set(gca,'xlim',xlim,'ylim',ylim);
  xlabel('training data')


  drawnow;
end;



N=size(data,1);% N=10;
perm=randperm(N);
data=data(perm,:);
targets=targets(perm);
countx=ones(1,N);
sumx=data(1:N,:)';

% function [countx,sumx]=clusters(countx,sumx,ncl);
% sequential clustering using Gini's measure 
% = vector quantization with constant scalar covariance
% provides only the cluster counts and sums
% which can then be used for classification purposes, using 
%     dim=size(sumx,1);
%     mu=sumc./countc(ones(dim,1),1:ncl)
% (perhaps after further coarsening by clusterh.m)
% countx(1,g)         number of vectors (or sum of weights) in group g
% sumx(:,g)           (weighted) sum of vectors in group g
% ncl                 number of final clusters 
disp('***********************************************************');
disp('***********************************************************');
disp('******************** before classification ****************');
disp('***********************************************************');
disp('***********************************************************');
disp([showtime,'   before classification']);
if sup>=0,
  % sequential clustering
  ncl=c;  % number of clusters
  [countc,sumc]=clusters(countx,sumx,ncl);
  dim=size(sumx,1);
  mu=sumc./countc(ones(dim,1),1:ncl)
  disp([showtime,' for clustering']);

  V=zeros(ncl,1);
  cpred=zeros(N,1);
  for l=1:N,
    for k=1:ncl,
      V(k)=sum((sumx(:,l)-mu(:,k)).^2);
    end;
    [Vmin,cpred(l,1)]=min(V);
  end;
end;
disp('***********************************************************');
disp('***********************************************************');
disp('******************** after classification *****************');
disp('***********************************************************');
disp('***********************************************************');

confusion=getfreq([targets,cpred])'
disp('rows=clusters found, columns=known classes');


if draw,
  subplot(1,2,2);cla;hold on;
  set(gca,'xlim',xlim,'ylim',ylim);
  for cl=1:ncl,
    ind=find(cpred==cl);
    plot(data(ind,1),data(ind,2),[type(cl),color(cl)]);
  end;
  xlabel('results of clustering')
end;

