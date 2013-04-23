%SYS10TOX convert from 10 based number to sys based
% sys and e have their least significant digit to the left
%[ e ] = sys10toX( index, sys, d )
% IN:
%   number(l,1)     l-th number.
%   sys             either: number, then conversion to [sys,..,sys] (d times)
%                   or:     [sys1,..,sysd]
%   d               desired dimension, can be omitted in case 2 above
% OUT:
%   e=[e1,..ed]     
%   fac, shift      index = e*fac-shift

% number has to be of type double! int types dont work (bug?)

function [ e, fac, shift ] = sys10toX( number, sys, d )

index = double(number) - 1;

n = length(index);

if length(sys) == 1,    
    sys = zeros(1,d) + sys;
else
    d = length(sys);
end

e = zeros(n,d);

for i=1:d
    e(:,i) = floor(mod(index,prod(sys(1:i))) / prod(sys(1:i-1))) + 1;
end    

if (nargout == 1), return; end

fac = ones(d,1);
shift = 0;

for i=1:d-1,
    fac(i+1) = sys(i)*fac(i);
    shift = shift + fac(i+1);
end