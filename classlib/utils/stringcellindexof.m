function [ idx ] = stringcellindexof( stringcellarray, string )
%STRINGCELLINDEXOF Summary of this function goes here
%   Detailed explanation goes here

for i=1:length(stringcellarray),
    if strcmp(stringcellarray{i},string),
        idx = i;
        return;
    end
end

idx = -1;

end

