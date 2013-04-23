% SUFF_NUM2POT Generate numerical potential information from a sufficient statistic.
%
% pot = suff_num2pot( suff, reg )
% 
% INPUT: 
%   suff    Sufficient statistic.
%   reg     (optional) regularization weight. 0 for no regularization. If >
%           0, total (over all classes) covariance estimation is added to 
%           each class covariance estimation with a weight of reg.
%
% OUTPUT:
%   pot(cl)  potential information struct for class cl with fields:
%       .freq   number of points.
%   	.mean   = sum(x)/.freq, mean estimation for this class.
%   	.cov    = sum(x*x')/.freq - mean*mean', covariance estimation for this class.
%   	.L      = (chol(.cov)')^-1, inverse lower cholesky factor of covariance matrix.
%   	.ent    = log(det(.cov)), class entropy.
%
% See also: suff_cat2pot, suff_pot.

function [ pot ] = suff_num2pot( suff, reg )

moment = suff.moment;
empty = [];

dim = size(moment,1) - 1;
ncl = size(moment,3);

if reg > 0,
    momtot = sum(moment,3);    
    potreg.cov = momtot(1:end-1,1:end-1)/momtot(end,end);
end

for c = 1:ncl,
    pot(c).freq = moment(end,end,c);
    
    if (pot(c).freq == 0)        
        pot(c).mean = zeros(dim,1);
        pot(c).cov = eye(dim);
        pot(c).L = eye(dim);
        pot(c).ent = inf;
        empty = [empty, c];
        continue;
    end
    pot(c).mean = moment(end,1:end-1,c)' / pot(c).freq; % nx1
    pot(c).cov = moment(1:end-1,1:end-1,c);

    if reg, % regularization w/ potential -> add it
        pot(c).cov = pot(c).cov + reg * potreg.cov;
    end

    % form covariance
    pot(c).cov = pot(c).cov / pot(c).freq - pot(c).mean*pot(c).mean';

    % regularize
    for i=1:dim,
        pot(c).cov(i,i)=pot(c).cov(i,i)*(1+sqrt(eps));
        if pot(c).cov(i,i)==0, pot(c).cov(i,i)=sqrt(realmin); end;
    end

    pot(c).L = inv(chol(pot(c).cov)');
    pot(c).ent = log(det(pot(c).cov));
end

if ~isempty(empty) % my favorite :)
    dprintf('suff_num2pot.m: there were empty classes:');
    empty
end