

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clusters.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [weight,sumdata,D,mu]=clusters(weight,sumdata,ncl);
% sequential clustering using Gini's measure 
% = vector quantization with constant scalar covariance
%
% weight(1,g)         number of vectors (or sum of weights) in group g
% sumdata(:,g)        (weighted) sum of vectors in group g
% ncl                 number of final clusters 
% D(f,g)              cost of merging final clusters f and g
% mu(:,g)             weighted mean of vectors in final cluster g   
%
function [weight,sumdata,D,mu]=clusters(weight,sumdata,ncl);

[dim,nall]=size(sumdata);
if size(weight,1)>1 | size(weight,2)~=nall,
  csize=size(weight),ssize=size(sumdata)
  error('input sizes do not match');
end;
if sum(weight<=0)>0, 
  error('counts must be positive');
end;
if ncl>nall,
  ncl,nall
  error('fewer groups than clusters');
end;
if ncl==nall,
  % nothing to be done
  return;
end;


% initialize cost matrix of size one higher 
% with first ncl+1 groups
ncl=ncl+1;

ni=1./weight(1:ncl);          %  1/n_k      (vector)
Ni=ni(ones(dim,1),:);         %  1/n_k      (matrix)
mu=sumdata(:,1:ncl).*Ni;         %  s_k/n_k    (vector)
G=mu'*mu;                     %  mu_k^Tmu_l
g=diag(G);                    %  ||mu_k||^2 (vector)
on=ones(1,ncl);
D=g(:,on);                    %  ||mu_k||^2 (matrix)
Ni=ni(on,:);                  %  1/n_k      (matrix)
D=D+D';D=(D-2*G)./(Ni+Ni');   %  ||mu_k-mu_l||^2n_kn_l/(n_k+n_l)
D=D+diag(inf*on);             % set diagonal to infinity

for next=nall:-1:ncl,
  if rem(next,1000)==0, disp(next), end;
  % find merge information
  [clist,klist]=min(D,[],1);
  [c,k]=min(clist);
  l=klist(min(find(clist==c)));

  % merge l into k
  weight(k)=weight(k)+weight(l);
  sumdata(:,k)=sumdata(:,k)+sumdata(:,l);

  % update cost matrix for row k
  ni(k)=1/weight(k);
  mu(:,k)=sumdata(:,k)/weight(k);
  ss=zeros(1,ncl);
  for s=1:dim,
    ss=ss+(mu(s,:)-mu(s,k)).^2;
  end;
  D(k,:)=ss./(ni+ni(k));
  D(:,k)=D(k,:)';
  D(k,k)=inf;

  % fill the free cluster 
  if l<next,
    weight(l)=weight(next);
    sumdata(:,l)=sumdata(:,next);
  end;
  if next==ncl, 
    ncl=ncl-1;
    weight=weight(1:ncl);
    sumdata=sumdata(:,1:ncl);
    % sort results by class size
    [weight,perm]=sort(-weight);weight=-weight;
    sumdata=sumdata(:,perm);
    D=D(perm,perm);
    % calculate cluster means
    dim=size(sumdata,1);
    mu=sumdata./weight(ones(dim,1),1:ncl);
    return; 
  end;

  % update cost matrix for row l
  k=l;
  ni(k)=1/weight(k);
  mu(:,k)=sumdata(:,k)/weight(k);
  ss=zeros(1,ncl);
  for s=1:dim,
    ss=ss+(mu(s,:)-mu(s,k)).^2;
  end;
  D(k,:)=ss./(ni+ni(k));
  D(:,k)=D(k,:)';
  D(k,k)=inf;

end;






