

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% recluster_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test recluster.m


% random data in 2 dimensions
N=200000;
ncl=32;
dim=2;
maxit=10;     % number of refinement iterations 
             % 1,4,10,inf are interesting
ncnew=15;
ndraw=1000;
markersize=4;

new=1;

if new,
  sumdata=fix(999.99999*rand(dim,N));
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


% weight(g)             sum of weights in group g
% sumdata(:,g)          weighted sum of vectors in group g
% g2cl(1,g)             initial cluster to which group g is assigned
%                       (if too many clusters, use instead cluster.m)
% percent               percentage of maximal cost allowed in merge
%
% g2cl(1,g)             final cluster to which group g is assigned
% mu(:,g)               weighted mean of vectors in final cluster g 
% cl2c(cl,nc)           class of cluster cl for a total of nc clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);clf;axis off;drawnow
figure(1);clf;axis off;drawnow
[g2cl,mu,cl2c]=recluster(weight,sumdata,g2cl,maxit,ncnew);

% input('scatter?>');
figure(2);clf;axis off;drawnow
scatterplot(data,g2cl,[],[],ndraw, markersize);
