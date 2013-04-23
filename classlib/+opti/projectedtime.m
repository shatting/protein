function [ etime ] = projectedtime( naa, sourcenaa, sourcetime )
%PROJECTEDTIME Summary of this function goes here
%   Detailed explanation goes here

etime = (2*naa-5).^3/((2*sourcenaa-5)^3/sourcetime);

end

