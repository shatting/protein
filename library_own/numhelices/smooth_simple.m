function [ y ] = smooth_simple( x, niters, kernelfn )
%SMOOTH_SIMPLE filter: smooth by neighbours (excl. point
% y(:,i) are smoothed function values after iteration i-1 (y(:,1) == x))
% kernelfn is the function used to determine value at
% y(i):=kernelfn(y(i-1),y(i+1))
% eg. simplesmooth(vals, 100, @(x,y)(x+y)/2) for mean or
% knum = simplesmooth(snum,100,@(x,y) max(x,y)) for maximum
% used for numsclasseslengthn.m

n=length(x);
y = zeros(n,niters+1);
y(:,1) = x;

for j=1:niters+1,
    tmp = y(:,j);
    for i=2:n-1,
    
        y(i,j+1) = kernelfn(tmp(i-1),tmp(i+1));
    
    end
end

end
