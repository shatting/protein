function [ count ] = countfrags( angle )
%COUNTFRAGS Summary of this function goes here
%   Detailed explanation goes here
count = 0;
for i=1:24,
    for j= 1:24,
        for k=1:24,
             if (~isempty(angle{i,j,k}))
                 count = count + 1;                 
             end
        end
    end
end