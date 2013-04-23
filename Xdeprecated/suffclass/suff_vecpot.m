%SUFF_VECPOT Calculate potential value of a (mixed) feature vector.
%
% [ v, cl, V ] = SUFF_VECPOT( x, pot, options )
% Calculate the (minimum) potential value(s) of a feature vector.
%
% INPUT:
%   x       the feature vector
%   pot     the potential (@see suff_pot.m)
%   options (optional)
%           .cat_ent_weight (default 1)
%           .num_ent_weight (default 1)
% OUTPUT:
%   [v, cl] = min(V) and V(c) = pot.cat(x,c) + pot.num(x,c) for all c
% 
% See also: suff_data2cl, suff_data2v.

function [ v, cl, V ] = suff_vecpot( x, pot, options )

if ~isfield(options,'cat_ent_weight'), options.cat_ent_weight = 1; end;
if ~isfield(options,'num_ent_weight'), options.num_ent_weight = 1; end;

catfeat = x(pot.cat);
numfeat = x(pot.num);
ncl = pot.ncl;

V = zeros(ncl,1);

if (isfield(pot,'cat'))
    if (isfield(pot.catpot,'lookup') && isfield(pot.catpot.lookup,'v'))
        ymax = pot.typ(pot.cat);
        xs = repmat(catfeat,ncl,1);
        iv = sysXto10([xs [1:ncl]'],[ymax ncl]);
        V = pot.catpot.lookup.v(iv);
    else
        V = suff_cat_vecpot(catfeat, pot.catpot, options.cat_ent_weight);
    end
end

if (isfield(pot,'num')),
    for c=1:ncl,      
        res=pot.numpot(c).L*(numfeat'-pot.numpot(c).mean);        
        V(c)= V(c) + res'*res - options.num_ent_weight*pot.numpot(c).ent;
    end    
end

[v, cl] = min(V);
