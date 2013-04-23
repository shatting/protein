function [ newt ] = oldt2newt( oldt )
%OLDT2NEWT Summary of this function goes here
%   Detailed explanation goes here

newt = atan2(oldt(:,2),oldt(:,1))/pi*100;

end
