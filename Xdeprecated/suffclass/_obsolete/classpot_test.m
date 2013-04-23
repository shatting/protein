%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% classpot_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test classpot.m, suffstat.m
% categorical/ordinal part not yet tested!
%


clear
clear mex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new=0;                  % new problem?
sup=1;                  % supervised (classes given)?
cat=0;                  % number of categories per feature 
                        % (0 = numerical features)
d=2;                    % dimension of data
c=5;                    % number of classes
np=[0,2,1,2,1]*50;      % number of data points per class
var=[1,1 2 2 1]*10;     % std of classes (>10 needed for rounding!)
fac=[1,1,1,1,1]*120;    % approx box size (fac/var=12 is interesting)

alg = 2;

% set parameters for classifier
opt=4;                  % only available option
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cat==1,error('cat must be 0 or >1'); end;
draw=( d==2 & c<=25 & ~cat );
type='o+x.*o+x.*o+x.*o+x.*o+x.*';
color='rrrrrbbbbbgggggmmmmmcccccyyyyy'; % for many classes
color='rbkmcyrbkmcyrbkmcyrbkmcyrbkmcy'; % for few classes

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
  oldalg = alg;
  load test_debug
  load test_debug0
  alg = oldalg;
end;
disp([showtime,' for data generation']);

mag=1;
data=mag*data;


if draw,
  rad=2.5; % factor for ellipse drawing; best in [2,3]
  % plot data 
  figure; clf; 
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

if sup,
  % supervised classification with numerical data
  ncl=c;
  typ=[typ,c];
  data=[data,targets];
  in=[1:d];out=[d+1];
else
  % unsupervised classification with numerical data
  % data already complete
  in=[1:d];out=in;
end;



% function suff=suffstat(ncl,typ);
%          initialize sufficient statistic 
%          for data of arbitrary type (>1 feature) and ncl classes
% function suff=suffstat(suff,data,N,class);
%          increment sufficient statistic using data(1:N,:) 
%          with assigned classes class(1:N)
% function [pot,potmin]=classpot(suff,opt);
%          create potential from sufficient statistic 
% function [cpred,prob]=classpot(pot,data,N,pr);
%          classify data(1:N,:) using the potential 
% function [dectree,pot,potmin]=classdec(typ,data,N,in,out,options);
%          create decision tree data of type typ
%          using data(1:N,in) to predict data(1:N,out) 
%          The potential pot (see classpot.m for details) allows one
%          to translate class information into data information
% function cpred=classdec(dectree,data,N);
%          classify incomplete data(1:N,:) using the decision tree
%          [this should not be used without preliminary pruning
%           but preferably postprocess with a pure potential method]

disp('***********************************************************');
disp('***********************************************************');
disp('******************** before classification ****************');
disp('***********************************************************');
disp('***********************************************************');
disp([showtime,'   before classification']);
if sup==1,
  % supervised classification with potential
  N=size(data,1);
  pr=np(1:c); % prior
  suff=suffstat(ncl,typ(in));
  suff=suffstat(suff,data(:,1:d),N,targets(:,1));
  disp([showtime,' for sufficient statistic']);
  if (alg == 1)
      [pot,potmin]=classpot(suff,opt);  
      disp([showtime,' for potential construction']);
      [cpred,prob]=classpot(pot,data,N,pr);
  else
      pot = suff_pot(suff);
      disp([showtime,' for potential construction']);
      cpred = suff_data2v(data,pot,pr);
  end
  disp([showtime,' for testing on training data']);
end;
disp('***********************************************************');
disp('***********************************************************');
disp('******************** after classification *****************');
disp('***********************************************************');
disp('***********************************************************');


if ~draw,
  freq=getfreq([targets,cpred'])
  disp('rows=classes, columns=leafs');
  disp('no figure')
end;

if draw,
  subplot(2,2,2);cla;hold on;
  set(gca,'xlim',xlim,'ylim',ylim);
  for cl=1:c,
    ind=(cpred==cl);
    plot(data(ind,1),data(ind,2),[type(cl),color(cl)]);
    if cat==0,
      % estimate ellipse
      if (alg == 1)
          Cest=inv(pot.R(:,:,cl)+pot.R(:,:,cl)')/pot.beta(cl);
          muest=pot.mu(:,cl);
      else
        Cest = inv(pot.numpot(cl).L * pot.numpot(cl).L');
        muest = pot.numpot(cl).mean;
      end
      ellipse(muest,rad,Cest,[color(cl),'-']);
    end;
  end;
  xlabel('results on training data')
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
    if (alg == 1)
        [cpred,prob]=classpot(pot,td,Ntest,pr);
    else
        cpred = suff_data2v(td,pot,pr);
    end
    for ic=c:-1:1,
      ind=(cpred==ic);
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
