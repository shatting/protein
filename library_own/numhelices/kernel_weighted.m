function [ y ] = kernel_weighted( x )
%KERNEL_WEIGHTED 
%   Detailed explanation goes here

weights = [1/length(x):1/length(x):1];
s = sort(x)';
y = sum(s.*weights)/length(x);

end
