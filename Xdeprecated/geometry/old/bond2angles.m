function [ angles ] = bond2angles( bond )
%bond2angles( bond ) returns dot products of consecutive bond vectors

n = length(bond);

angles = zeros(n-1,1);

for i=1:length(angles),
    angles(i) = bond(i,:)*bond(i+1,:)';
end