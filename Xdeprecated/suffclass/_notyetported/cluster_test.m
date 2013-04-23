

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cluster_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test cluster.m


% random data in 2 dimensions
N=10000;
ncl=50;
dim=2;

setting.nrand=3;         % number of random clusterings
setting.ncluster=10;      % number of initial clusters wanted
setting.ncnew=6;         % number of clusters wanted after reclustering
setting.maxit=10;         % maximal number of refinement iterations
                         % 1,4,10,inf are interesting
setting.draw=[1 2];      % index set for drawing scatterplots
setting.ndraw=500;       % maximal number of points drawn per cluster
setting.markersize=3;    % markersize for scatterplot

new=1;      % new data? (0=repeat old data)

if new,
  sumdata=fix(999.99999*rand(2,N));
  ind=find(abs(sum(sumdata,1)-1000)>100);
  N=length(ind);
  sumdata=sumdata(:,ind);
  data=sumdata/1000;
  weight=1000*ones(1,N);
  
  g2cl=1+rem(randperm(N),ncl);
  g2clold=g2cl;
  Nold=N;
else
  N=Nold;
  g2cl=g2clold;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g2cl,mu,cl2c]=cluster(weight,sumdata,setting);

