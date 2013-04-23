%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% clusterh_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test clusterh.m
%


clear
clear mex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new=1;                  % new problem?
cat=10;                  % number of categories per feature 
                        % (0 = numerical features)
d=2;                    % dimension of data
c=4;                    % number of classes
np=[2,1,2,1]*50;        % number of data points per class
var=[1 2 2 1]*10;       % std of classes (>10 needed for rounding!)
fac=[1,1,1,1]*200;      % approx box size 
                        % (fac/var=12 is interesting, >12 easier)


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
  figure(1); clf; 
  subplot(2,2,1);hold on;
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
countx=ones(1,N);
sumx=data(1:N,:)';

% function [merge,cost,sep]=clusterh(countx,sumx);
% hierarchical clustering using Gini's measure 
% = vector quantization with constant scalar covariance
%
% countx(1,g)         number of vectors (or sum of weights) in group g
% sumx(:,g)           (weighted) sum of vectors in group g
% merge(:,i)=[l;k]    merge in step i group l into group k<l 
% cost(1,i)           cost of merging in step i
% sep                 separating index set at the root (=1)     
disp('***********************************************************');
disp('***********************************************************');
disp('******************** before classification ****************');
disp('***********************************************************');
disp('***********************************************************');
disp([showtime,'   before classification']);
if sup>=0,
  % hierarchical clustering
  [merge,cost,sep]=clusterh(countx,sumx);
  disp([showtime,' for clustering']);

  ncl=c+1; % show one leaf more than classes
  cpred=zeros(N,1);
  cpred(1)=1;
  for i=1:ncl-1,
    l=merge(1,N-i); 
    cpred(l)=i+1;
  end;
  for i=N-ncl:-1:1,
    l=merge(1,i); 
    k=merge(2,i);
    cpred(l)=cpred(k);
  end;
end;
disp('***********************************************************');
disp('***********************************************************');
disp('******************** after classification *****************');
disp('***********************************************************');
disp('***********************************************************');

endcost=cost(end-20:end)/max(cost)
subplot(2,2,3);cla
semilogy(cost)

confusion=getfreq([targets,cpred])'
disp('rows=clusters found, columns=known classes');


if draw,
  subplot(2,2,2);cla;hold on;
  set(gca,'xlim',xlim,'ylim',ylim);
  for cl=1:ncl,
    ind=find(cpred==cl);
    plot(data(ind,1),data(ind,2),[type(cl),color(cl)]);
  end;
  xlabel('results of clustering')
end;

