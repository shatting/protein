

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mom2pot.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [meanx,Rx,beta,Vmin,sigx]=mom2pot(countx,sumx,Mx,reg);
% function [meanx,Rx,beta,Vmin,sigx]=mom2pot(suff,reg);
%
% create potential coefficients from moment information
% generated elsewhere from data vectors with a class assignment
%
% function mom2pot(prt);  % for debug
% repeat calculation with previous data and printlevel prt
%  
% The potential has the form
%    V(x,t)=beta(t)*(x-meanx(:,t))'*Rx(:,:,t)*(x-meanx(:,t))+Vmin(t)
% with integral Rx, and is related to the Gaussian densities by
%    p(x|t)=exp(-V(x,t)).
% For prediction with prior class probabilities pr(t),
% minimize V(t|x)=V(x,t)-log(pr(t))
%
% The covariance matrices can be obtained by
%    C=inv(Rx(:,:,t)+Rx(:,:,t)')/beta(t);
% the potential mean is
%    Vmean(t)=Vmin(t)+nd/2
% where nd=size(sumx,1) is the dimension of the data vectors.
%
% Regularization ensures that R+R' is positive definite 
% 
% countx(1,t)     number of vectors of class t
% sumx(:,t)       sum  of vectors of class t
% Mx(:,:,t)       2nd moment matrix of vectors of class t
% suff            sufficient statistic containing countx,sumx,Mx
% reg             number of dummy entries, must be >1 (default 2)
%
% meanx(:,t)      mean for class t (rounded to integer)
%                 or empty (if data set is too small; 
%                           then Rx contains the error message)
% sigx(:,t)       standard deviations for class t (optional)
% Rx(:,:,t)       upper triangular quadratic term for class t
% beta(1,t)       potential scaling factor
% Vmin(1,t)       potential minimum
%                 
function [meanx,Rx,beta,Vmin,sigx]=mom2pot(countx,sumx,Mx,reg)

% extract moment information
if nargin==2, reg=sumx; end;
if nargin<=2,
    suff=countx;
    countx=squeeze(suff.moment(end,end,:));
    sumx=squeeze(suff.moment(1:end-1,end,:));
    Mx=squeeze(suff.moment(1:end-1,1:end-1,:));
end;

if nargin<4 && nargin ~=2, reg=2; end;
if reg<=1, reg=1+eps; end;

% create total mean and total covariance matrix 
[nd,nc]=size(sumx);
if nd==0, 
  Rx='no numerical features present';
  meanx=[];beta=[];Vmin=[];sigx=[];
end;

nn=0;
mu=zeros(nd,1);
C=zeros(nd,nd);
for t=1:nc,
  nt=countx(t);
  if nt>0,
    nn=nn+nt;
    mu=mu+sumx(:,t);
    C=C+Mx(:,:,t)-sumx(:,t)*sumx(:,t)'/nt;
  end;
end;
if nn<=nd, 
  Rx='too few data for meaningful covariance matrix';
  meanx=[];beta=[];Vmin=[];sigx=[];
end;
mu=mu/nn;
C=C/(nn-1);
% take care of null rows
normC=norm(C(:),inf);
nullrows=find(diag(C)<=0)';
dC=diag(C);dC(nullrows)=1+dC(nullrows);
C=C+diag(dC/nn); 
M=C+mu*mu'; % moment contribution per dummy entry

% create regularized class means, covariance matrices
% and initial potential
% V=fac*[(x-mu)^TC^(-1)(x-mu)+log det(2pi C)]
meanx=zeros(nd,nc);
sigx=zeros(nd,nc);
Rx=zeros(nd,nd,nc);
Vmin=zeros(1,nc);
beta=zeros(1,nc);

for t=1:nc,
  % add reg dummy entries
  nt=countx(t)+reg;
  st=sumx(:,t)+reg*mu;
  Ct=Mx(:,:,t)+reg*M;
  
  % round mean
  %%mut=round(st/nt);
  mut=st/nt;
  
  Ct=(Ct-mut*st'+(nt*mut-st)*mut')/(nt-1);
  normCt=norm(Ct(:)); 
  if normCt>0,
    beta(t)=1/(2*sqrt(nt*nd)*normCt);
    H=inv(Ct)/beta(t);
    badH=(sum(isinf(H(:)))>0);
  else
    badH=1;
  end;
  if badH, 
    error('singular Ct');
  end;

  % round quadratic term 
  % unrounded: D=diag(diag(H));R=triu(H-D);D=0.5*D;R=D+R;H=R+R';
  D=diag(diag(H));
  %R=round(triu(H-D)); % off-diagonal part
  R=triu(H-D); % off-diagonal part
  lam=0.5*max(real(eig(H-R-R',D))); % approx 1/2
  %D=fix(lam*D)+eye(nd);
  D=lam*D+eye(nd);
  R=D+R;
  H=R+R';

  % constant term
  Vmin(t)=0.5*sum(log((2*pi/beta(t))./abs(eig(H))));

  % assign to table
  Rx(:,:,t)=R;
  meanx(:,t)=mut; 

  if nargout==5,
    Cte=inv(H)/beta(t);
    sigx(:,t)=(1+fix(10*sqrt(diag(Cte))))/10;
  end;

  if 0, % debug mode
    Ct
    Cte=inv(H)/beta(t)
    lC=0.5*log(det(2*pi*Ct))
    lCe=0.5*log(det(2*pi*Cte))
    Vmint=Vmin(t)
    % check quality of rounding 
    % relerr should be around 0.2 but never close to 1
    relerr=(trace(R*Ct)*beta(t)-0.5*nd)/sqrt(2*nd/nt)
  end;

end;

