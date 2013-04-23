function [ y ] = shift( x, v )
%SHIFT shifts all vectors x(i,:) or x(:,i) by the vector v
% if v is a column vector, x is assumed to be of the second form

if (size(v,1)<size(v,2))
    y = x + repmat(v,size(x,1),1);
else
    y = x + repmat(v,1,size(x,2));
end