function [ oldt ] = newt2oldt( newt )
%OLDT2NEWT Summary of this function goes here
%   Detailed explanation goes here

% newt \in [-100,100]
% oldt \in [-1,1]x[-1,1]

oldt(:,1) = cos(double(newt)/100*pi);
oldt(:,2) = sin(double(newt)/100*pi);

end
