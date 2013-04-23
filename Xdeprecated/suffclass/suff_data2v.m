% SUFF_DATA2V Vectorized class prediction and potential evaluation.
%
% [ cl, vmin, V, prob ] = SUFF_DATA2V( data, pot, pr, options )
% Calculate class assignment, potential values and conditional
% probabilities for data under potential pot. This is the vectorized
% version of suff_data2cl.
%
% INPUT:
% pot           potential (suff_pot.m)
% data(l,:)     l-th data vector
% pr            (optional) input priors. if not given, then uniform assumed
% options       .ent_weight (default 1)
%               .use_entropy (default false)    if true, pr is assumed to
%                                               be entropy values instead
%                                               of probabilities.
%               .handleequalpotentials          false, "warn", "error", "unclassified" or
%                                               "randomize" (TODO:
%                                               partially implemented)
%
% OUTPUT: 
% cl(l,1)       predicted class for lth data vetor
% vmin(l,1)     minimum potential value
% V(l,cl)       potential value of lth vector under potential cl
% prob(l,cl)    probability that lth vector is in class cl
%
% See also: suff_pot, suff_data2cl.

function [ cl, vmin, V, prob, Vdiff ] = suff_data2v( data, pot, pr, options )

if nargin<4,
    options = struct;
end
if ~isfield(options,'ent_weight'), options.ent_weight = 1; end
if ~isfield(options,'use_entropy'), options.use_entropy = 0; end
if ~isfield(options,'handleequalpotentials'), options.handleequalpotentials= 0; end

nnum = length(pot.num);
ncat = length(pot.cat);

if nargin>=3 && ~isempty(pr),
  if length(pr)~=pot.ncl,
    pot.ncl,size(pr) 
    error('prior length must equal number of classes');
  end;  
end;
for i=1:ncat,
    ii=pot.cat(i);
    if min(data(:,ii))<=0, 
      error('categorical variables must be positive')
    else
      [dmax,ind]=max(data(:,ii));
      if dmax>abs(pot.typ(ii)), 
        ind,ii,dmax
        error('data(ind,ii) exceeds number of levels')
      end;
    end;
end;

N=size(data,1);
on=ones(N,1);


for cl=pot.ncl:-1:1, 
  %disp(['get potential values for class ',num2str(cl)]);
  if nargin==2 || sum(pr) == 0,
    % use uniform prior
    potcl=zeros(N,1);
  else
    % use input prior or entropy
    if (options.use_entropy),
        potcl=zeros(N,1) + options.ent_weight*pr(cl);
    else
        %TODO: either consensus classification or figure out a way to add
        %cat and num potentials without adding prior 2 times
        %potcl=zeros(N,1) - options.ent_weight*log(pr(cl)+realmin);
        potcl = zeros(N,1) + options.ent_weight*pot.catpot.ent(cl);
    end
  end;
  
  if nnum>0,
    row=pot.numpot(cl).mean';
    diff=data(:,pot.num)-row(on,:);
    diffR=diff*pot.numpot(cl).L';
    for k=1:nnum,
      potcl=potcl+diffR(:,k).*diffR(:,k);
    end;
%    potcl=pot.beta(cl)*potcl;
%    potcl=potcl+pot.shift(1,cl);
  end;
  
  if ncat>0,
    potcl=potcl';
    for i=1:ncat,
      potcl=potcl+pot.catpot.pot(i,data(:,pot.cat(i)),cl);
    end;
    potcl=potcl';
  end;
  V(:,cl)=potcl;
  
end;
[vmin,cl]=min(V,[],2); %TODO: equal potentials: randomize class assignment?

if options.handleequalpotentials
    Vdiff = sort(abs(V - repmat(vmin,1,pot.ncl)),2); % sort potential differences; then the first column will always be zero
    Vdiff = min(Vdiff(:,2:end),[],2);
    Vequal = Vdiff < 1e-10;
    if any(Vequal),
        warning(sprintf('%i (%d%%) of data points had two or more classes with the same potential value.',sum(Vequal),sum(Vequal)/length(Vequal)*100));
    end
end

if nargout<4, return; end;
% get probability table
for cli=pot.ncl:-1:1, 
  epot(:,cli)=exp(vmin-V(:,cli)); 
end;
epotsum = sum(epot,2);
prob=epot./epotsum(:,ones(pot.ncl,1));
