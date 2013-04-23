
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% covclass.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pot,cl,confus]=covclass(data,target,weight);
% partially supervised classification of continuous data 
% target labels are known classes if positive, tentative if negative, 
% unknown if 0
%
% The predicted classes may have different labels than initially.
% Predict original target labels by label=argmin(cost*confus(:,cl))
% where cost(s,t) is the cost of deciding for label s given target t
%
% data(l,:)             l-th data vector
% target(l,1)           l-th class number 
%                       fixed if positive, tentative if negative, 
%                       unknown if 0
% weight(l,1)           weight of data vector
% 
% pot                   potential defining classes by minimizing
%                       V_c(x) = ||L(x-mean)||^2+ent 
%                       over the classes with  V_c(x)<=Vmax_c 
%   .freq               class frequency
%   .mean               class mean 
%   .cov                class covariance
%   .ent                class entropy
%   .L                  class transform
%   .Vmax               maximal reliable potential value
% cl(l,1)               final class assignment
% confus(t,c)           number of labels +t in class c
%                       (but if all t<=0, of labels -t)
%
function [pot,cl,confus]=covclass(data,target,weight);

if nargin==3, error('weights not yet implemented'); end;
tim=toc;

% initial class assignment
[nd,dim]=size(data)    % nd items of length dim
if nd==0, error('empty data set'); end;
data=double(data);
if min(target)<0, 
   cl=abs(double(target));
else
   cl=double(target);
end;
if any(cl==0), 
  ncl=max(cl)+1;   % number of classes
  cl(cl==0)=ncl; 
else
  ncl=max(cl);   % number of classes
end; 
freq=zeros(1,ncl);     
for c=1:ncl,
  freq(c)=sum(cl==c);
end;

maxitbad=10;
level(1:10)=inf;
iteration=0;
while 1,
  iteration=iteration+1

  % sort classes by frequency
  [nfsort,perm]=sort(freq);
  invp(perm)=[ncl:-1:1]';
  cl=invp(cl);
  % merge tiny classes 
  large=(sum(freq.^2>=nd));
  if large<ncl, 
    ncl=1+sum(freq.^2>=nd);
    cl(cl>ncl)=ncl;
  end;

  % find sufficient statistic
  pot.freq=zeros(1,ncl);      % class frequency
  pot.mean=zeros(dim,ncl);    % class mean mu
  pot.cov=cell(1,ncl);        % class covariance C
  pot.L=cell(1,ncl);          % class transform, L^-1 L^-T = C
  pot.ent=zeros(1,ncl);       % class entropy
  pot.Vmax=zeros(1,ncl);      % class potential maximum
  for c=1:ncl,
    pot.cov{c}=zeros(dim,dim);
    pot.L{c}=zeros(dim,dim);
  end;
  for c=1:ncl,
    ind=(cl==c);
    cdata=data(ind,:);
    pot.freq(c)=size(cdata,1);
    pot.mean(:,c)=sum(cdata,1);
    pot.cov{c}=cdata'*cdata;
  end;

  % get total covariance model
  tfreq=sum(pot.freq);
  tmean=sum(pot.mean,2)/tfreq;
  tcov=pot.cov{1};
  for c=2:ncl,
    tcov=tcov+pot.cov{c};
  end;
  tcov=tcov/tfreq-tmean*tmean';
  % regularization
  for i=1:dim,
    tcov(i,i)=tcov(i,i)*(1+sqrt(eps));
    if tcov(i,i)==0, tcov(i,i)=sqrt(realmin); end;
  end;

  % get class covariance matrices
  for c=ncl:-1:1,
    if c==ncl,
      % add mean as sample point
      pot.freq(c)=pot.freq(c)+1;
      pot.mean(:,c)=pot.mean(:,c)+tmean;
      pot.cov{c}=pot.cov{c}+tmean*tmean';
    end
    % regularize second moments by adding one total covariance
    cov= pot.cov{c}+tcov;
    % form mean and covariance
    pot.mean(:,c)=pot.mean(:,c)/pot.freq(c);
    cov=cov/pot.freq(c)-pot.mean(:,c)*pot.mean(:,c)';
    for i=1:dim,
      cov(i,i)=cov(i,i)*(1+sqrt(eps));
    end;
    pot.cov{c}=cov;
    L=chol(cov)';
    pot.L{c} = inv(L);
    pot.ent(c) = log(det(L));
  end;

  % get new class predictions
  V=zeros(1,ncl);
  Vd=zeros(nd,1);
  nclold=ncl;
  clold=cl;cl=0*cl;
  tim1=toc;
  for l=1:nd,
    togo=nd+1-l;
    NN=1000;
    if l>NN & rem(togo,NN)==0,
      timetogo=showtime((toc-tim1)*togo/(l-1));
      disp([num2str(togo),' items to come in approx. ',timetogo]);
    end;
    x=data(l,:)';
    for c=1:ncl,
      res=pot.L{c}*(x-pot.mean(:,c));
      V(c)= res'*res+pot.ent(c);
    end; 
    [Vd(l),cl(l)]=min(V);
  end;
  freq=zeros(1,ncl);     
  for c=1:ncl,
     freq(c)=sum(cl==c);
  end;
  freq(ncl)=freq(ncl)+1; % regularization

  % stopping test
  ndif=sum(clold~=cl);
  disp(' ');
    old_new_frequencies=[pot.freq;freq]
  disp('note one dummy item for the mean');
  level
  disp([num2str(ndif),'/',num2str(nd),' reassignments made']);
  showtime(toc-tim)
  disp(' ');
  if ndif==0 | ndif>=max(level), break; end;
  level=[level(2:10),ndif];
  input('next>');
end;
  

% confus(c,cl)          number of labels +c in class cl; 

nct=max(target);
if nct==0, 
  target=-target;
  nct=max(target);
end;
nct=double(nct);
confus=zeros(nct,ncl);
for t=1:nct,
  clist=cl(target==t);
  for c=1:ncl,
    confus(t,c)=sum(clist==c);
  end;
end;

disp('pot.Vmax not yet implemented');




