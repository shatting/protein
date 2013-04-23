function [ cl, m ] = applymerge( cl, merge )
%APPLYMERGE Summary of this function goes here
%   Detailed explanation goes here

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
m(setdiff(m,cls)) = 0; %     1    2    3    0    0    4
m = m';
cl = m(cl);

end
