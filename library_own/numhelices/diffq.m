function [ x ] = diffq( x )
%DIFFQ Summary of this function goes here
%   Detailed explanation goes here

x = (x(3:end)-x(2:end-1))./(x(2:end-1)-x(1:end-2));

end
