%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% classpot.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pot,potmin,freq]=classpot(suff,opt);
%          create potential from sufficient statistic 
%          rounding data to integers is assumed to be be harmless!! 
% function [cpred,prob,datapot]=classpot(pot,data,N,pr,options);
%          classify data(1:N,:) using the potential pot and priors pr  
%
% suff       sufficient statistic
%            .ncl      number of classes
%            .typ      specifies type of variables
%            .ncat     number of categorical variables
%            .nnum     number of numerical variables
%            .nord     number of ordinal variables
%            .num      list of numerical features
%            .cat      list of categorical features
%            .ord      sublist of ordinal features
%                      ord(cat) is the list of ordinal features
%            .moment   symmetric moment matrix for numerical features
%                      moment(:,:,cl)= sum [x 1]^T*[x 1] 
%                      over all x=data(l,num) with class(l)=cl
%                      moment(end,end,:) is the class size
%                                        even if only categorical data
%            .count    count for categorical features
%                      count(i,a,cl) = number of l with 
%                                      class(l)=cl, data(l,cat(i))=a
% opt        defines recipe for potential creation (default; 4)
%               opt = 1: C_i = C for all classes i
%               opt = 2: C_i = sigma(i)*I (I unity matrix)
%               opt = 3: C_i = diag(sigma(i,:))
%               opt = 4: general C_i
% pot        information defining the potential
%            V(data) = Vnum(data(num)) + sum_i ccat(data(cat(i)))
%            Vnum(x) = (x-mu)'*R*(x-mu) + shift
%            the estimated conditional probability is
%            prob(data|cl)=exp(-V(data)) 
%            the total distribution needs a class prior
%            .ncl      number of classes
%            .typ      specifies type of variables
%            .ncat     number of categorical variables
%            .nnum     number of numerical variables
%            .cat      list of categorical features
%            .num      list of numerical features
%            .ccat     coefficients for categorical potential
%                      ccat(i,a,cl) for class cl, data(cat(i))=a
%            .mu       center vectors for numerical potential 
%                      mu(:,cl) for class cl
%            .sig      standard deviations for numerical potential
%                      sig(:,cl) for class cl
%            .R        quadratic term for numerical potential
%                      R(:,:,cl) (upper triangular) for class cl
%            .beta     scale factors for numerical potential
%            .shift    shifts for numerical potential
%                      shift(1,cl) for class cl
% data(l,:)  l-th feature vector (a row)
% N          number of vectors to be classified
% pr(cl)     prior probability or frequency of class cl
%            (default: all-one)
% cpred(l,1) predicted class assignment of l-th feature vector
% prob(l,cl) probability that l-th feature vector is in class cl
% potmin(cl) potential minimum of class cl 
% freq(cl)   regularized relative frequency of class cl in training

function [cpred,prob,freq]=classpot(p,data,N,pr,options)

if (nargin < 5)
    options = struct;
end

if (~isfield(options,'cat_ent_weight')) options.cat_ent_weight = 1; end

if nargin<3 || isempty(data),
  % function [pot,potmin,freq]=classpot(suff,opt);
  % create potential from sufficient statistic 
  % p contains suff, data contains recipe, 
  % pot is returned in cpred, potmin in prob
  if nargin==1, opt=4;
  else          opt=data;
  end;
  reg=1+p.nnum;  % regularization factor; required reg>1
  countx=p.moment(end,end,:);
  countall=sum(countx);
  if countall==0, error('no data present'); end;
  if opt==4,
    if p.nnum>0,
      ind=1:size(p.moment,1)-1;
      sumx=p.moment(ind,end,:);
      Mx=p.moment(ind,ind,:);
      [meanx,Rx,beta,Vmin,sig]=mom2pot(countx,sumx,Mx,reg);
      if isempty(meanx),error(Rx); end;
    else
      meanx=[];Rx=[];beta=[];Vmin=[];sig=[];
    end;
  else
    error(['option ',num2str(opt),' not yet implemented']);
  end;

  if p.ncat>0,
    on=ones(1,size(p.count,2),1);
    maxcat = max(p.typ(p.cat));
    regterm=(1./p.typ(p.cat)')*reg(on);
    % set positions correctly
    shorts = find(p.typ(p.cat)< maxcat);
    for i=shorts,
        regterm(i,p.typ(p.cat(i))+1:maxcat) = 0;
    end
    
    for cl=1:p.ncl,
      countcl=p.count(:,:,cl)+regterm;
      %logterm=zeros(size(countcl))-inf;
      %ind=find(countcl>0);
      %logterm(ind)=-log(countcl(ind)./countsum(ind));
      %ccat(:,:,cl)=log(freq(cl))+logterm;
      countsum=sum(countcl,2);
      ccat(:,:,cl)=-log(countcl./countsum(:,on));
    end;
  else
     ccat=[];
  end;

  % assign potential (in cpred) 
  p=struct('ncl',p.ncl,'typ',p.typ,'ncat',p.ncat,'nnum',p.nnum,...
           'cat',p.cat,'num',p.num,'ccat',ccat,...
           'mu',meanx,'sig',sig,'R',Rx,'beta',beta,'shift',Vmin);
  cpred=p;

  if nargout<2, return; end;
  
  % assign potmin (in prob)
  prob=zeros(1,p.ncl);
  for cl=1:p.ncl,
    if p.nnum>0, prob(cl)=p.shift(cl);
    else         prob(cl)=0;
    end;
    if p.ncat>0,
      prob(cl)=prob(cl)+min(min(p.ccat(:,:,cl)));
    end;
  end;

  % assign regularized relative frequency
  countx=countx+reg;
  countall=countall+double(p.ncl)*reg;
  freq=squeeze(countx)/countall;

  return;
end;




% function [cpred,prob]=classpot(p,data,N,pr);
% classify data(1:N,:) using the potential p and the prior pr
% V(data) = Vnum(data(num)) + sum_i ccat(data(cat(i)))
% Vnum(x) = beta*(x-mu)'*R*(x-mu) + shift

if nargin==4,
  if length(pr)~=p.ncl,
    p.ncl,size(pr) 
    error('prior length must equal number of classes');
  end;
  for i=1:p.ncat,
    ii=p.cat(i);
    if min(data(:,ii))<=0, 
      error('categorical variables must be positive')
    else
      [dmax,ind]=max(data(:,ii));
      if dmax>abs(p.typ(ii)), 
        ind,ii,dmax
        error('data(ind,ii) exceeds number of levels')
      end;
    end;
  end;
end;

data(N+1:end,:)=[];
on=ones(N,1);
for cl=p.ncl:-1:1, 
  %disp(['get potential values for class ',num2str(cl)]);
  if nargin==3 || sum(pr) == 0,
    % use uniform prior
    potcl=zeros(N,1);
  else
    % use input prior
    potcl=zeros(N,1)-options.cat_ent_weight*log(pr(cl)+realmin);
  end;
  
  if p.nnum>0,
    row=p.mu(:,cl)';
    diff=data(:,p.num)-row(on,:);
    diffR=diff*p.R(:,:,cl);
    for k=1:p.nnum,
      potcl=potcl+diffR(:,k).*diff(:,k);
    end;
    potcl=p.beta(cl)*potcl;
    potcl=potcl+p.shift(1,cl);
  end;
  
  if p.ncat>0,
    potcl=potcl';
    for i=1:p.ncat,
      potcl=potcl+p.ccat(i,data(:,p.cat(i)),cl);
    end;
    potcl=potcl';
  end;
  pot(:,cl)=potcl;
  
end;
[potmin,cpred]=min(pot');
if nargout==1, return; end;

freq = pot;

% get probability table
potmin=min(pot,[],2);
for cl=p.ncl:-1:1, 
  epot(:,cl)=exp(potmin-pot(:,cl)); 
end;
epotsum=sum(epot')';
prob=epot./epotsum(:,ones(p.ncl,1));

% debug information
return;
round(100*[pot/max(pot(:)),potmin/max(pot(:)),cpred',prob])
max(pot(:))
p,p.R
error('stop')
