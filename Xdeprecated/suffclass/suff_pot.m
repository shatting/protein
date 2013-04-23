% SUFF_POT Generate potential information from sufficient statistic.
%
% pot = SUFF_POT ( suff, options )
% create potential from sufficient statistic.
%
% INPUT:
%   suff        sufficient statistic (see suffstat)
%   options     (optional) 
%               .cat_cond_reg = realmin
%               .cat_prior_reg = 0.1
%               .cat_ent_weight = 0.5
%               .cat_no_lookup = 1
%               .reg_weight = 1
%
% OUTPUT:
%  pot      potential struct with fields
%               .numpot     see suff_num2pot.m    
%               .catpot     see suff_cat2pot.m
%
% See also: suffstat, suff_num2pot, suff_cat2pot.

function [ pot ] = suff_pot( suff, options )

if (nargin<2), options = struct; end;
if ~isfield(options,'cat_cond_reg'), options.cat_cond_reg = realmin; end;
if ~isfield(options,'cat_prior_reg'), options.cat_prior_reg = 0.1; end;
if ~isfield(options,'cat_ent_weight'), options.cat_ent_weight = 0.5; end;
if ~isfield(options,'cat_no_lookup'), options.cat_no_lookup = 1; end;
if ~isfield(options,'reg_weight'), options.reg_weight = 1; end;

% numeric portion
if (suff.nnum > 0),
    pot.numpot = suff_num2pot(suff, options.reg_weight);
end

% categorical portion
if (suff.ncat > 0),    
        
    pot.catpot = suff_cat2pot(suff, options.cat_prior_reg, options.cat_cond_reg);    
    
    if (~options.cat_no_lookup && suff.nnum == 0)
        try
            pot.catpot.lookup = suff_catpot2lookup(suff.typ(suff.cat), pot.catpot, 1, options.cat_prior_reg);
        catch
            dprintf('not using cat lookup, feature space is too big (%i items)',prod(suff.typ(suff.cat)));
        end
    end
end

pot.ncl = suff.ncl;
pot.num = suff.num;
pot.cat = suff.cat;
pot.typ = suff.typ;
  