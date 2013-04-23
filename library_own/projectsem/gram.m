function [ g ] = gram( bond )
%GETANGLES Summary of this function goes here
%   Detailed explanation goes here


n=length(bond);
angles = zeros(n-1,1);

for i=1:length(bond)-1,

  angles(i)=bond(i,:)*bond(i+1,:)';

end

