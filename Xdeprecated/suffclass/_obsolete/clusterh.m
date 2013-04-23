

% handling of tiny classes for split not yet correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clusterh.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [merge,cost,D]=clusterh(countx,sumx);
% hierarchical clustering using Gini's measure 
% = vector quantization with constant scalar covariance
%
% function [merge,cost,sep]=clusterh(countx,sumx,splitlim);
% split a set into two by hierarchical clustering
% 
% countx(1,g)         number of vectors (or sum of weights) in group g
% sumx(:,g)           (weighted) sum of vectors in group g
% splitlim            minimal cluster weight for split (default 0)
%
% merge(:,i)=[l;k]    merge in step i group l into group k<l 
% cost(1,i)           cost of merging in step i
% D(f,g)              cost of merging initial groups f and g
% sep                 index list for a separating cluster 
%                     (next to root if splitlim<=1) 
%                     weight and complementary weight >=splitlim
%                     if this is impossible, sep=[]
%
function [merge,cost,sep]=clusterh(countx,sumx,splitlim);

[dim,ng0]=size(sumx);
if size(countx,1)>1 | size(countx,2)~=ng0,
  csize=size(countx),ssize=size(sumx)
  error('input sizes do not match');
end;

% handle zero weights
ind=find(countx>0);
ng=length(ind);
if ng<=1,
  % no split
  merge=[];cost=[];sep=[];
  return;
end;
countx(countx<=0)=NaN+countx(countx<=0); 


merge=zeros(2,ng-1);
cost=zeros(1,ng-1);

% initialize cost matrix

ni=1./countx;                 %  1/n_k      (vector)
Ni=ni(ones(dim,1),:);         %  1/n_k      (matrix)
mu=sumx.*Ni;                  %  s_k/n_k    (vector)
G=mu'*mu;                     %  mu_k^Tmu_l
g=diag(G);                    %  ||mu_k||^2 (vector)
on=ones(1,ng0);
D=g(:,on);                    %  ||mu_k||^2 (matrix)
Ni=ni(on,:);                  %  1/n_k      (matrix)
D=D+D';D=(D-2*G)./(Ni+Ni');   %  ||mu_k-mu_l||^2n_kn_l/(n_k+n_l)
D=D+diag(inf*on);             % set diagonal to infinity

if nargin<3,
  % store distance matrix in third output argument
  sep=D;
end;

for i=1:ng-1,
  % find merge information
  [clist,klist]=min(D(ind,ind),[],1);
  [c,k]=min(clist);
  l=min(find(clist==c));
  k=ind(k);
  l=ind(klist(l));
  if k<l, merge(1,i)=l;merge(2,i)=k;
  else    merge(1,i)=l;merge(2,i)=k; % interchange would cause trouble?
  end;
  cost(i)=c;
  % update cost matrix 
  % perhaps slightly cheaper by restricting to ind?
  countx(k)=countx(k)+countx(l);
  ni(k)=1/countx(k);
  sumx(:,k)=sumx(:,k)+sumx(:,l);
  mu(:,k)=sumx(:,k)/countx(k);
  ss=zeros(1,ng0);
  for s=1:dim,
    ss=ss+(mu(s,:)-mu(s,k)).^2;
  end;
  D(k,:)=ss./(ni+ni(k));
  D(:,k)=D(k,:)';
  D(k,k)=inf;
  ind(ind==l)=[];
end;

if nargin<3, return; end;



% find separating index set
if nargin<3, splitlim=0; end;
countall=countx(ind); % total sum of weights
splitlim2=countall-splitlim; 
kok=0;
lok=0;
splitok=0;
ind=zeros(1,ng0);
for i=ng-1:-1:1,
  l=merge(1,i);
  k=merge(2,i);
  countx(k)=countx(k)-countx(l);
  if ~splitok,
    kok=(countx(k)>=splitlim & countx(k)<=splitlim2);
    lok=(countx(l)>=splitlim & countx(l)<=splitlim2);
    splitok=(kok|lok);
    if kok, ind(k)=1;kok=countx(k); end;
    if lok, ind(l)=2;lok=countx(l); end;
  else,
    ind(l)=ind(k);
  end;
end;
if max(kok,lok)==0, sep=[];           wt=0;
elseif kok>=lok,    sep=find(ind==1); wt=kok;
else                sep=find(ind==2); wt=lok;
end;

% split=[wt,countall-wt]

