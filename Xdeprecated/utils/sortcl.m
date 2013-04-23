function [ cl ] = sortcl( cl )
%SORTCL sort classes by frequency

[nfsort,perm]=sort(freqs(cl));
invp(perm)=[max(perm):-1:1]';
cl=invp(cl)';

end
