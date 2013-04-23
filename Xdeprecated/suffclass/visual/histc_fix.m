%HISTC_FIX Slightly altered version of built-in histc.
%Adds last bin count to second-to-last and removes last. 
%(histc makes extra count for last)
%
% See also: histc.

function [ n ] = histc_fix( x, edges )

n = histc(x,edges);
n(end-1) = n(end-1) + n(end);
n(end) = [];

end
