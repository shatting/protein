%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% classdec_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test classdec.m
%


clear
clear mex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new=1;                  % new problem?
sup=2;                  % 0: unsupervised (clustering)
                        % 1: supervised (classes given)
                        % 2: supervised (numerical values given)
cat=20;                  % number of categories per feature 
                        % (0 = numerical features)
d=2;                    % dimension of data
c=4;                    % number of classes
np=[2,1,2,1]*200;        % number of data points per class
var=[1 2 2 1]*10;       % std of classes (>10 needed for rounding!)
fac=[1,1,1,1]*400;      % approx box size 
                        % (fac/var=12 is interesting, >12 easier)


% set parameters for classifier
opt=4;                  % only available option
prt=0;                  % printlevel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cat==1,error('cat must be 0 or >1'); end;
draw=( d==2 & c<=35 & ~cat );
type ='o+x.*o+x.*o+x.*o+x.*o+x.*o+x.*o+x.*';
color='rrrrrbbbbbkkkkkmmmmmcccccyyyyyggggg';
color='rbkmcygrbkmcygrbkmcygrbkmcygrbkmcyg'; % preferable

save test_debug0 new sup opt draw

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
  data=round(data);

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

if sup==1,
  % supervised classification with categorical target
  typ=[typ,c];
  data=[data,targets];
  in=[1:d];out=[d+1];
elseif sup==2,
  % supervised classification with numerical target
  typ=[typ,0];
  data=[data,targets];
  in=[1:d];out=[d+1];
else
  % unsupervised classification (target ignored)
  % data already complete
  in=[1:d];out=in;
end;
N=size(data,1);
disp([showtime,'   before classification; reset to zero']);



% function [dectree,pot,pmin,cpred]=classdec(typ,data,N,in,out,options);
%          create decision tree data of type typ
%          using data(1:N,in) to predict data(1:N,out) 
%          The potential pot (see classpot.m for details) allows one
%          to translate class information into data information
% function cpred=classdec(dectree,data,N);
%          classify incomplete data(1:N,:) using the decision tree
%          [this should not be used without preliminary pruning
%           but preferably postprocess with a pure potential method]
% options  options influencing the tree construction
%          options(1)=opt defines recipe in potential training 
%	     opt = 1: C_i = C for all classes i
%	     opt = 2: C_i = sigma(i)*I (I unity matrix)
%	     opt = 3: C_i = diag(sigma(i,:))
%	     opt = 4: general C_i
%          options(2)=clsiz: minimal class size
%          options(3)=prc: percentage of data declared outliers
%          options(5)=cbysort: handling of categorical questions
%                     1: sort by min potential
%                     0: use simple clustering
options(1)=4;
options(2)=max(d+1,fix(N/12));
options(3)=0;
options(4)=prt;

fprintf(1,'\n\n\n\n\n\n');
options(5)=0;

disp('***********************************************************');
disp('******************** before classification ****************');
disp('***********************************************************');
tic
if sup>=0,
  % classification with decision tree
  [dectree,pot,pmin,cpred]=classdec(typ,data,N,in,out,options);
  disp([showtime,' for training dectree']);
  cpred2=classdec(dectree,data,N);
  disp([showtime,' for testing on training data']);
end;
disp('***********************************************************');
disp('******************** after classification *****************');
disp('***********************************************************');


if ~draw,
  freq=getfreq([targets,cpred]);
  disp('rows=classes, columns=leafs');
  freq2=sortrows(freq')'
  % disp('same table with sorted columns')
  %freq3=getfreq([targets,cpred2])
  %disp('table recomputed from decision tree');
  disp('no figure')
end;

if draw,
  subplot(2,2,2);cla;hold on;
  set(gca,'xlim',xlim,'ylim',ylim);
  ncl=dectree(1);
  for cl=1:ncl,
    ind=find(cpred==cl);
    plot(data(ind,1),data(ind,2),[type(cl),color(cl)]);
  end;
  xlabel('results on training data (from trainer)')

  subplot(2,2,3);cla;hold on;
  set(gca,'xlim',xlim,'ylim',ylim);
  ncl=dectree(1);
  for cl=1:ncl,
    ind=find(cpred2==cl);
    plot(data(ind,1),data(ind,2),[type(cl),color(cl)]);
  end;
  xlabel('results on training data (from classifier)')
  drawnow;
end;

if draw,
  for pic=0:1,
    subplot(2,2,3+pic);cla;hold on;
    sc=3*pic; % scale of test picture
    xlim=xlim+(xlim(2)-xlim(1))*[-sc,sc];
    ylim=ylim+(ylim(2)-ylim(1))*[-sc,sc];
    set(gca,'xlim',xlim,'ylim',ylim);
    % create test data on a grid
    ngrid=137;
    Ntest=(1+ngrid)^2;
    xgrid=ones(ngrid+1,1)*[0:ngrid]/ngrid;
    ygrid=xgrid';
    td(1:Ntest,1)=xlim(1)+(xlim(2)-xlim(1))*xgrid(:);
    td(1:Ntest,2)=ylim(1)+(ylim(2)-ylim(1))*ygrid(:);
    if cat>0,
      for i=1:2,
        td(:,i)=td(:,i)-min(td(:,i));
        td(:,i)=round(0.501+(cat-0.002)*(td(:,i)/max(td(:,i))));
      end;
    end;
    cpred2=classdec(dectree,td,Ntest);
    for ic=ncl:-1:1,
      ind=(cpred2==ic);
      h=plot(td(ind,1),td(ind,2),['.',color(ic)]);
      set(h,'markersize',9);
    end;
    xlabel('results on test data')
    drawnow;
  end;
end;

showtime(toc)


return;

% old stuff
suff=suffstat(1,typ(in));
ind=find(targets==1);N1=size(ind,1)
suff=suffstat(suff,data(ind,1:d),N1,targets(ind,1));
[pot1,potmin1]=classpot(suff,opt);

suff=suffstat(1,typ(in));
ind=find(targets==2);N1=size(ind,1)
suff=suffstat(suff,data(ind,1:d),N1,targets(ind,1)-1);
[pot2,potmin2]=classpot(suff,opt);
pmin=[potmin;potmin1,potmin2]/mag^2
pmu=[pot.mu;pot1.mu,pot2.mu]/mag
pR=[pot.R(:,:,1),pot.R(:,:,2);pot1.R,pot2.R]
