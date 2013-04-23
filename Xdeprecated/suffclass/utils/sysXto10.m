function [ n, fac, shift ] = sysXto10( v, sys )
% [ n, fac, shift ] = sysXto10( v, sys )
% v is n 1xp numbers in 1xp system sys
% inverse of v = sys10toX(n,[sys])
%error if only one class at the end

v = double(v); sys = double(sys);
fac = ones(1,length(sys));
for i=1:length(sys)-1,
    fac(i+1) = prod(sys(1:i));    
end

shift = sum(fac(2:end));

n = (fac*v' - shift)';