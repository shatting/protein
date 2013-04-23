function [ x ] = removezeros( x, kernelfn )
%REMOVEZEROS filter: removes zeros
% set all zeros z to kernelfn(z-,z+) where z- = value of next non-zero to
% the left and z+ value of next non-zero to the right
% used in numsclasseslengthn.m 

if (nargin < 2), kernelfn = @(x,y) (x+y)/2; end

n=length(x);

i=min(find(x~=0)) + 1;
while i<n,
    if (x(i) == 0)
        lowi = max(find(x(1:i)~=0));
        hii = min(find(x(i:end)~=0)) + i - 1;
        x(lowi+1:hii-1) = kernelfn(x(lowi),x(hii));
                
        i = hii + 1;
    else
        i = i+1;
    end    
end


end
