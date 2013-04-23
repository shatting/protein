% Multinomial Naive Bayes (independent occurrence assumption)
% 
% function [w] = nbMulti(x,y,alpha,twoclass)
% x - examples [n,d]
% y - integer labels (1..l) [n,1]
% alpha - smoothing constant
% twoclass - if non-zero, assume one-versus-all binary task with
%            'twoclass' being the positive class [scalar]
% w - learned weight vector [d,l] ([d,1] if twoclass is non-zero)
% 
% Written by Jason Rennie, January 2004
% Last modified: Tue May 24 17:53:55 2005

function [w] = nbMulti(x,y,alpha,twoclass)
  fn = mfilename;
  if (nargin < 3 | nargin > 4)
    error('%s: incorrect number of arguments',fn);
  end
  if nargin == 3
    twoclass = 0;
  end
  [n,d] = size(x);
  len = sum(x,2);
  if twoclass
    pidx = find(y==twoclass);
    nidx = find(y~=twoclass);
    ptot = sum(len(pidx));
    ntot = sum(len(nidx));
    pcnt = full(sum(x(pidx,:),1))';
    ncnt = full(sum(x(nidx,:),1))';
    w = log((pcnt+alpha)./(ptot+d*alpha)) - log((ncnt+alpha)./(ntot+d*alpha));
  else
    for i=1:max(y)
      idx = find(y==i);
      tot = sum(len(idx));
      cnt = full(sum(x(idx,:),1))';
      w(:,i) = log((cnt+alpha)./(tot+d*alpha));
    end
  end


