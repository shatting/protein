function [t] = cell2structarray( data, naa )
%cell2structarray extract a random struct array of naa chains from
%DS.mat-like cell array

if (nargin==2)
    nchains = min(naa,length(data));
    perm = randperm(length(data));
else
    nchains = length(data);
    perm = [1:nchains];
end

j = 1;

perm = sort(perm(1:nchains));

for i=perm,
   t(j).res = data{i}.res;
   t(j).name = data{i}.name;
   t(j).seq = data{i}.seq;
   t(j).bond = data{i}.bond;
   
   j = j + 1;
end