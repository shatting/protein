function [ res ] = findprotswithlength( data, length )
%findprotswithlength( data, length ) find all proteins by length 

n = size(data,2);
res = [];
for i=1:n,
    if (size(data{i}.bond,1) + 1 == length)
        res=[res, i];
	end
end