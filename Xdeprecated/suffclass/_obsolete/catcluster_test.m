

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% catcluster_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test cat*.m categorical clustering

clear
clear mex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new=1;                  % new problem?
cat=4;                  % number of categories per feature 
d=2;                    % dimension of data
c=4;                    % number of classes
np=[2,1,2,1]*5;         % number of data points per class
var=[1 2 2 1]*10;       % std of classes (>10 needed for rounding!)
fac=[1,1,1,1]*200;      % approx box size 
                        % (fac/var=12 is interesting, >12 easier)


% set parameters for classifier
prt=1;                  % printlevel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cat<=1,error('cat must be >1'); end;

save test_debug0 new 

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


[data,targets]

disp('not yet programmed');

return;
confusion=getfreq([targets,cpred])'
disp('rows=clusters found, columns=known classes');



