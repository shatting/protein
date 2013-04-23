%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cluster.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [g2cl,mu,cl2c]=cluster(weight,sumdata,setting);
% clustering using Gini's measure 
% = vector quantization with constant scalar covariance
% calls the other clustering routines
%
% nrand random sequential clusterings 
% are followed by a consensus clustering
% relevant dendrograms are plotted
%
% weight(g)             number of vectors (or sum of weights) in group g
% sumdata(:,g)          (weighted) sum of vectors in group g
%                       sumdata(g,:) is allowed if sumdata is not square
% setting               parameter setting
%   .nrand                 number of random clusterings
%   .ncluster              number of initial clusters wanted
%   .ncnew                 number of clusters wanted after reclustering
%   .maxit                 maximal number of refinement iterations
%   .draw                  index set for drawing scatterplots
%   .ndraw                 maximal number of points drawn per cluster
%   .markersize            markersize for scatterplot
%
% g2cl(1,g)             final cluster to which group g is assigned
% mu(:,g)               weighted mean of vectors in final cluster g   
% cl2c(cl,nc)           class of cluster cl for a total of nc clusters
%                       (if fewer clusters are desired)
%
function [g2cl,mu,cl2c]=cluster(weight,sumdata,setting);

nrand=setting.nrand;
ncluster=setting.ncluster;
maxit=setting.maxit;
ncnew=setting.ncnew;
draw=setting.draw;
ndraw=setting.ndraw;
markersize=setting.markersize;

% transpose if necessary

N=length(weight);
[dim,N1]=size(sumdata);
if N~=N1,
  if dim==N,
    sumdata=sumdata';
    dim=N1;
  else
    N,dim,size(weight)
    error('sizes do not fit');
  end;
end;
if N~=size(weight,2), 
  weight=weight'; 
end;
mutot=sum(sumdata,2)/sum(weight);
on=ones(dim,1);
mudata=sumdata./weight(on,:);

figure(1);clf;axis off;drawnow
nrow=fix(sqrt(nrand+1));
ncol=fix((nrand+1)/nrow);
if nrow*ncol<=nrand,ncol=ncol+1; end;
mulist=zeros(dim,ncluster,nrand);
clweight=zeros(nrand,ncluster);
clabel=zeros(N,nrand);
for cas=1:nrand,cas

  % random sequential clustering
  perm=randperm(length(weight));
  countx=weight(perm);sumx=sumdata(:,perm);
  [countc,sumc,D]=clusters(countx,sumx,ncluster);
  mu=sumc./countc(ones(dim,1),:);

  % iterative refinement of classes
  [clabel(:,cas),mu,cweight]=clusteri(weight,mudata,mu,maxit);
  mulist(1:size(mu,1),1:size(mu,2),cas)=mu;
  clweight(cas,1:size(mu,2))=cweight;
  sumx=mu.*cweight(ones(dim,1),1:ncluster);
  [merge,cost,D]=clusterh(cweight,sumx);

  % plot results
  figure(1);subplot(nrow,ncol,cas)
  dendro(cweight,merge,cost,D);
  figure(cas+1);clf;axis off;drawnow
  scatterplot(mudata(draw,:),clabel(:,cas),[],[],ndraw, markersize);
  drawnow
  %input('cont>');
end;
np1=min(20,dim);
np2=min(13,ncluster);
sizes=clweight(:,1:np2)
mutot_mu=[mutot(1:np1,1,ones(1,nrand)),mulist(1:np1,1:np2,:)]



disp('reclustering');
 

% find coarsest refinement of all clusterings
cmax=max(clabel);
cltot=clabel(:,1)-1;
for cas=2:nrand,
  cltot=cltot*cmax(cas-1)+clabel(:,cas)-1;
end;
cltot=cltot+1;
g2cl=0*cltot;
cnew=0;
clist=rsort(cltot)';
clist=clist(clist>0);
for cc=clist
  ind=find(cltot==cc);
  nc=length(ind);
  if nc>0,
    cnew=cnew+1;
    g2cl(ind)=cnew+g2cl(ind);
  end;
end;

% get clustering at prescribed percent level
figure(1);subplot(nrow,ncol,nrand+1);axis off;
[g2cl,mu,cl2c]=recluster(weight,sumdata,g2cl,maxit,ncnew);

figure(nrand+2);clf;axis off;drawnow
scatterplot(mudata(draw,:),g2cl,[],[],ndraw, markersize);


