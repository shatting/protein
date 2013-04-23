% SUFF_CAT_VECPOT Calculate potential value of a categorical feature vector.
% This function is very slow.
%
% v = suff_cat_vecpot( x, catpot, entfactor )
%
% INPUT: 
%   x           1xd categorical feature vector.
%   catpot      categorical potential (see suff_cat2pot)
%   entfactor   factor to apply to potential entropies.
%
% OUTPUT:
%   v           nclx1 class potential values.
%
% See also: suff_vecpot, suff_cat2pot.

function v = suff_cat_vecpot( el, pot, entfactor )
   
if (isstruct(entfactor)), entfactor = entfactor.cat_ent_weight; end; % entfactor = options

d = size(pot.pot,1);
ncl = size(pot.pot,3);

% compute class potentials
if 0, % this is faster, but wrong
   p = pot.pot(:,el,:);
   p = reshape(p,d^2,ncl);   
   v = sum(p(1:d:d^2));       
else
   v = zeros(1,ncl);   
   for i=1:d,
        v = v + reshape(pot.pot(i,el(i),:),1,ncl);
   end
end

% add entropy
v = v + entfactor * pot.ent;

