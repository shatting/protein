%APPLYMERGE Apply merging information to a list of class numbers.
%
%[ clmerged, m ] = applymerge( cl, merge )
% INPUT:
%   cl(l,1)     l-th class
%   merge       in ith step, merge merge(1,i) into merge(2,i)
%
% OUTPUT: 
%   clmerged    = m(cl)
%   m           merge information after all steps
%
% See also: suff_clusterh.

function [ cl, m ] = applymerge( cl, merge )

if (isempty(merge)) return; end

mi = size(merge,2);

for i=1:mi,
    mg = merge(:,i);
    cl(cl==mg(1)) = mg(2);
end

% remove empty classes < maxclass
cls = unique(cl);
m = 1:max(cls);
m(cls) = find(cls);
m(setdiff(m,cls)) = 0;
m = m';
cl = m(cl);

end
