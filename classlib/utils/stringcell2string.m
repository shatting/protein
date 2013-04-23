function [ s ] = stringcell2string( sc, delimiter )
%STRINGCELLPRINT Summary of this function goes here
%   Detailed explanation goes here

if (nargin < 2)
    delimiter = ', ';
end

s = sc{1};
for i=2:length(sc),
    s = [s,delimiter,sc{i}];
end

end

