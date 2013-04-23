function [ t1, t2, t3, t4 ] = sys10to20( index )
%SYS10TO20 Summary of this function goes here
%   Detailed explanation goes here

index = index - 1;

t1 = floor(index / 20^3) + 1;
t2 = floor(mod(index,20^3) / 20^2) + 1;
t3 = floor(mod(index,20^2) / 20) + 1;
t4 = mod(index,20) + 1;