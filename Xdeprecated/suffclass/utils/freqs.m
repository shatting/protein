function [ freq ] = freqs( classes, ncl )
%freqs( classes ) class frequencies. can handle class 0 and empty classes.

if nargin < 2,
    ncl = max(classes);
end   

freq = zeros(ncl,1);
for i=1:length(freq),
    freq(i) = sum(classes == i);
end