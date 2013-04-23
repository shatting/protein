function [ index ] = sys20to10(t1,t2,t3,t4)
%SYS20TO10 Summary of this function goes here
%   Detailed explanation goes here

if (nargin == 1)
    t2 = t1(2);
    t3 = t1(3);
    t4 = t1(4);
    t1 = t1(1);
end

index = (t1-1)*20^3 + (t2-1)*20^2 + (t3-1)*20 + t4;