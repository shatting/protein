% SUFF_CATPOT2LOOKUP Compute lookup table for a categorical potential. 
% Done by iteration over all elements in feature space (slow).
%
% lookup = SUFF_CATPOT2LOOKUP( ymax, pot, saveclasspots, entfactor )
%
% lookup.v(y1,..,yd,c) := sum_i(pot.catpot(i,yi,c)) + entfactor*pot.ent(c)
% lookup.g(y1,..,yd)   := argmin_c [ sum_i(pot.catpot(i,yi,c)) + entfactor*pot.ent(c) ]
%
% Usage of lookup:
% clpred = lookup.cl(data*lookup.fac-lookup.shift);
%   where data(l,:) is l-th feature vector
%   and clpred(l) is predicted class of l-th data vector
%
% INPUT:
%   ymax(i)         number of cats in position i      
%   pot.pot(i,a,cl) potential for y_i=a in class cl
%   pot.ent(cl)     entropy of class cl
%   saveclasspots   boolean, all class potential values saved if true (default 0)
%   entfactor       (default 1)
%
% OUTPUT:
%   lookup  .cl(y1,..,yd)       class of y
%           .v(y1,..,yd,cl)     potential value of y in class cl
%           .fac  
%           .shift
%
% See also: suff_catpots.

function [ lookup ] = suff_catpot2lookup( ymax, pot, saveclasspots, entfactor )

if nargin < 3, saveclasspots = 0; end;
if nargin < 4, entfactor = 1; end;

ncl = size(pot.pot,3);
d = size(pot.pot,1);

% initialize lookups; can throw error if too big
lookup.cl = zeros(ymax);
if saveclasspots,
    try    
        lookup.v = zeros([ymax,ncl]); % could throw error even if .cl doesnt
        saveclasspots = 1;
    catch
        warning('categorical potential values not saved.');
    end
end

% get fac and shift
[el,lookup.fac,lookup.shift] = sys10toX(1,ymax);

% run through all possible feature combinations
% TODO: could be done faster&more elegantly i think
casmax = prod(ymax(1:end));
rep = fix(casmax/10);
for cas = 1:casmax,
   if (rem(cas,rep) == 0), dprintf('lookup %.2f%%',cas/casmax*100); end   
   
   % get element
   el = sys10toX(cas,ymax);
   
   % get class potentials
   v = suff_cat_vecpot(el, pot, entfactor);
   
   % save class potentials for element
    if saveclasspots
       els = repmat(el,ncl,1);
       iv = sysXto10([els [1:ncl]'],[ymax ncl]);
       lookup.v(iv) = v;
    end
    
   % get minimizing class
   [vs,lookup.cl(cas)] = min(v);
end