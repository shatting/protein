function [ idx ] = findchain( data, chainname )
%FINDCHAIN finds all chains with names starting with chainname

idx = [];
l = length(chainname);

for i=1:length(data),
    if strcmp(data{i}.name(1:l),chainname),
        idx = [idx i];
    end
end


end
