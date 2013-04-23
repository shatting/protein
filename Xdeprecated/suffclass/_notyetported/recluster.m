

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% recluster.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [g2cl,mu,cl2c]=recluster(weight,sumdata,g2cl,maxit,ncnew);
% reclustering using Gini's measure 
% = vector quantization with constant scalar covariance
% calls the other clustering routines
%
% weight(g)             sum of weights in group g
% sumdata(:,g)          weighted sum of vectors in group g
% g2cl(1,g)             initial cluster to which group g is assigned
%                       (if too many clusters, use insted cluster.m)
% maxit                 maximal number of refinement iterations
% ncnew                 merge to ncnew clusters (or fewer)
%
% g2cl(1,g)             final cluster to which group g is assigned
% mu(:,g)               weighted mean of vectors in final cluster g 
% cl2c(cl,nc)           class of cluster cl for a total of nc clusters
%
function [g2cl,mu,cl2c]=recluster(weight,sumdata,g2cl,maxit,ncnew);

cla;drawnow;
dim=size(sumdata,1);
mu=zeros(dim,0);
ng=0;
for c=1:max(g2cl),
  ind=find(g2cl==c);
  sumnew=sum(sumdata(:,ind),2);
  wtnew=sum(weight(ind));
  if wtnew>0,
    ng=ng+1;
    mu(:,ng)=sumnew/wtnew;
  end;
end;
on=ones(dim,1);
mudata=sumdata./weight(on,:);
[g2cl,mu,clweight]=clusteri(weight,mudata,mu,maxit);
ncl=length(clweight)
sumx=mu.*clweight(on,:);
[merge,cost,D]=clusterh(clweight,sumx);

% merge to retain ncnew clusters
ncl=max(g2cl);
ccurrent=[1:ncl]';
for i=1:ncl-ncnew,
  % merge in step i group l into group k<l 
  l=merge(1,i);
  k=merge(2,i);
  ind=(ccurrent==l);
  ccurrent(ind)=k+0*ccurrent(ind);
end;
% label new clusters as 1:nc
nc=0;
cl2c=zeros(ncl,1);
for c=1:ncl,
  ind=find(ccurrent==c);
  if ~isempty(ind),
    nc=nc+1;
    cl2c(ind)=nc+0*cl2c(ind);
  end;
end;
disp([num2str(ncl),' clusters reduced to ',num2str(nc)]);

% find new cluster means
sumnew=zeros(dim,nc);
wtnew=zeros(1,nc);
for cl=1:ncl;
  ind=find(g2cl==cl);
  c=cl2c(cl);
  sumnew(:,c)=sumnew(:,c)+sum(sumdata(:,ind),2);
  wtnew(:,c)=wtnew(:,c)+sum(weight(ind));
end;
mu=sumnew./wtnew(on,:);


% recluster with new means
[g2cl,mu,clweight]=clusteri(weight,mudata,mu,maxit);
on=ones(dim,1);
sumx=mu.*clweight(on,:);
[merge,cost,D]=clusterh(clweight,sumx);
cl2c=dendro(clweight,merge,cost,D);  

np1=min(20,dim);
ncl=size(mu,2);
np2=min(13,ncl);
ncl
sizes=clweight
mutot=sum(sumdata,2)./sum(weight);
mutot_mu=[mutot(1:np1,1),mu(1:np1,1:np2)]

