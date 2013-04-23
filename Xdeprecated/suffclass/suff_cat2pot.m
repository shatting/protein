% SUFF_CAT2POT Generate categorical potential information from a sufficient statistic.
% catpot = SUFF_CAT2POT( suff, rpriors, rconditionals ) 
%
% catpot.pot(i,a,cl) = -log(count(i,a,cl)+rconditionals)+log(freq(c)+rconditionals*d);
% catpot.ent(cl)     = -log(freq(cl) + rpriors);
% (d = length(suff.cat))
%
% INPUT:
%   suff        sufficient statistic
%   rpriors     number of prior regularization points 
%   rcond       number of conditionals regularization points
%
% OUTPUT:
%   catpot      .pot(i,a,cl)        potential for a in position i in class cl
%               .ent(cl)            cl-th class entropy    
%               .entsum             entsum(:,:,cl) == ent(cl). for simple
%                                   summation with .pot.
%                                   exp(-(pot.pot+pot.entsum)) ~ suff.count (~
%                                   b/c of regularization)
%               .freq               1xncl class freqs
%               .prob(i,a,cl)       P(C=cl|A_i = a)
%               .condprob(cl,a_1,..,a_d)     naive conditional probability
%                                   P( cl | A = (a_1,..,a_d) )
%
% @smoothing:
% http://www.cs.cmu.edu/~tom/mlbook/NBayesLogReg.pdf    
%
% See also: suffstat.

function pot = suff_cat2pot( suff, rpriors, rconditionals )

if (nargin < 2)
    rpriors = 0;
end
if (nargin < 3)
    rconditionals = 0;
end
dprintf('regularization with %g class points and %g conditionals.', rpriors, rconditionals);

count = suff.count;
ncl = suff.ncl;
ymax = suff.typ(suff.cat);
d = length(ymax);
maxcat = max(ymax);
countdimmask = zeros(d,maxcat,ncl);
for i=1:d, countdimmask(i,1:ymax(i),:) = 1; end;

% counttot(i,a,c) = #category a is on place i, independent of c ~ P(a_i=a)
counttot = sum(count,3);
counttot = counttot(:,:,ones(ncl,1)); 

% countclass(i,a,c) = #datapoints in class c
countclass = repmat(sum(sum(count,1),2)/d,size(count,1),size(count,2));

% .freq, class frequencies
pot.freq = squeeze(sum(suff.count,2));
pot.freq = pot.freq(1,:);

% .pot
if (rconditionals == 0)
    ind = count>0;
    z = count==0;
    nzeroconds = sum(z(:));
    if (nzeroconds > 0)
        dprintf('warning: number of zero conditional probabilities: %i. should use conditionals regularization.', nzeroconds);
    end
    pot.pot = zeros(size(count)) + inf;
    pot.pot(ind) = -log(count(ind)) + log(countclass(ind));
else	
    dimensioncats = repmat(ymax',1,[maxcat ncl]);% strange error here with matlab2009b: second dimension is 127 long instead of 400
    %reason: ncl is int8, and matlab strangely converts the array to int8
    %and not to double (as it should i think)
    %fix: make ncl double in suffclass. no, make classes double    
    pot.pot = -log(count+rconditionals) + log(countclass.*dimensioncats+rconditionals*dimensioncats);
    pot.pot(~countdimmask) = inf;
end

% .ent
if (rpriors == 0) 
   pot.ent = zeros(1,ncl) + inf;
   indf = pot.freq>0;
   nzeropriors = sum(pot.freq==0);    
   if (nzeropriors > 0)
       dprintf('warning: number of empty classes: %i. should use priors regularization.', nzeropriors);
   end
   pot.ent(indf) = -log(pot.freq(indf));
else
   pot.ent = -log(pot.freq + rpriors);
end
    
% .entsum
pot.entsum = repmat(pot.ent',1,[d maxcat]);
pot.entsum = permute(pot.entsum,[2 3 1]);

% .prob
ind = counttot > 0;
pot.prob = zeros(d,max(ymax),ncl);
pot.prob(ind) = count(ind)./counttot(ind);
