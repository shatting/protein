%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cat2group.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function gpred=cat2group(lookup,data);
% function [gpred,confusion]=cat2group(lookup,data,target);
% function [gpred,confusion,g2t,mc]=cat2group(lookup,data,target,cost);
%
% classification of categorical data by table lookup
%
% data(l,:)             l-th data vector
% lookup                group dictionary; 
% lookup.g(y1,y2,...)   group label of item (y1,y2,...)
% lookup.fac            index factors
% lookup.shift          index shift
%                       dft(y1,y2,...)=dft(y*lookup.fac-lookup.shift)
% target(l,1)           l-th target (= class number to be predicted)
% cost(s,t)             cost of deciding for target s given target t
%                       must be nonnegative with zero diagonal
%                       (default: 1-eye, also used if cost=1)
%
% gpred(l,1)            predicted group for l-th data vector
% confusion(t,g)        number of targets t in group g; 
% g2t(1,g)=t            least costly class assignment t to group g
% mc(1,g)               misclassification cost for this assignment
% 
% total misclassification cost:    mctot=sum(mc);
% number of misclassified items:   nmisc=sum(g2t(gpred)~=target);
%
function [gpred,confusion,g2t,mc]=cat2group(lookup,data,target,cost);


% compute group assignment
[N,d]=size(data);
if isa(data,'double'),
  gpred=lookup.g(data*lookup.fac-lookup.shift);
else
  % convert to double in batches
  gpred=uint8(zeros(N,1));
  batchsize=100000;
  nn=fix(batchsize/d);
  nb=fix(N/nn)+1;
  for ind1=1:nn:N,
    ind=[ind1:min(N,ind1+nn-1)];
    gpred(ind)=lookup.g(double(data(ind,:))*lookup.fac-lookup.shift);
  end;
end;
if nargout==1, return; end;

% compute confusion matrix
confusion=getfreq([target,gpred]);
if nargout==2, return; end;

% compute target prediction table
if nargin<4,cost=1;end;
if length(cost)==1,
  [n,nn]=size(confusion);
  if n~=nn, error('default size requires group=class'); end;
  cost=1-eye(n);
end;
[mc,g2t]=min(cost*confusion);g2t=g2t(:)';

