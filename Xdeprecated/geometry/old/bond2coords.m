function [ coords ] = bond2coords( bond )
%bond2coords( bond )

bond = double(bond);

coords = zeros(size(bond,1)+1,3);

for i=2:size(bond,1)+1,
    coords(i,:) = coords(i-1,:) + bond(i-1,:);
end