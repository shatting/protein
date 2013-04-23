%SORTCL Sort a list of numbers by frequency.
% 1 will have the highest frequency.

function [ cl ] = sortcl( cl )

[nfsort,perm]=sort(freqs(cl));
invp(perm)=[max(perm):-1:1]';
cl=invp(cl)';

end
