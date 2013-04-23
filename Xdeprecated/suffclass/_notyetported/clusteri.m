

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clusteri.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [clabel,munew,cweight]=clusteri(weight,mudata,mu,maxit);
% iterative clustering using Gini's measure 
% = vector quantization with constant scalar covariance
%
% weight(1,g)         number of vectors (or sum of weights) in group g
% mudata(:,g)         (weighted) mean of vectors in group g
% mu(:,c)             center of cluster c
% maxit               maximal number of iterations (default; inf)
%
% clabel(1,g)         cluster label for group g
% munew(:,c)          center of groups labeled by cluster c
% cweight(1,c)        weight of cluster c
%
function [clabel,mu,cweight]=clusteri(weight,mudata,mu,maxit);

prt=1;

if nargin<4, maxit=inf; end;
[dim,N]=size(mudata);
[dimc,nc]=size(mu);
if dim~=dimc, 
  dim,dimc
  error('sizes do not match'); 
end;
if size(weight,2)~=N, 
  N,size(weight,2)
  error('sizes do not match'); 
end;

oldlabel=[];
while maxit>0,
  maxit=maxit-1;

  % get new labeling
  V=inf+zeros(1,N);
  clabel=zeros(N,1);
  for g=1:nc,
    if N>100000, disp(['group ',num2str(g)]); end;
    Vg=0;
    for k=1:dim,
      Vg=Vg+(mudata(k,:)-mu(k,g)).^2;
    end;
    ind=find(Vg<V);
    clabel(ind)=g+0*clabel(ind);
    V(ind)=Vg(ind);
  end;

  % update means
  on=ones(1,dim);
  for c=nc:-1:1,
    ind=find(clabel==c);
    weightc=sum(weight(ind));
    if weightc>0,
      mu(:,c)=sum(weight(on,ind).*mudata(:,ind),2)/weightc;
    end;
  end;

  % check confusion
  if ~isempty(oldlabel),
    confusion=getfreq([oldlabel,clabel])';
    nexchanged=sum(confusion(:))-sum(diag(confusion));
  else
    nexchanged=inf;
  end;
  if prt>1,
    if ~isempty(oldlabel),
      confusion
      disp('rows=clusters found, columns=old clusters');
      sizes_old_new=[sum(confusion,1);sum(confusion,2)']
    else
      cweight=[];
      for c=nc:-1:1,cweight(1,c)=sum(weight(clabel==c));end;
      cweight
    end;
  end;
  if nexchanged<=0, break; end;
  oldlabel=clabel;
  if prt & nexchanged<inf, 
    disp([,'exchanged ',num2str(nexchanged)]), 
  end;
  if prt>1, input('next iteration>'); end;
end;

cweight=[];
for c=nc:-1:1,cweight(1,c)=sum(weight(clabel==c));end;
% sort by increasing size
[cweight,perm]=sort(-cweight);cweight=-cweight;
nc=max(find(cweight>0));
cweight=cweight(1:nc);
perm=perm(1:nc);
mu=mu(:,perm);
newclass(perm)=[1:nc]';
clabel=newclass(clabel);
clabel=clabel(:);

