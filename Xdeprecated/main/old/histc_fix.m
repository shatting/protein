function [ n ] = histc_fix( x, edges )
%HISTC_FIX returns counts between edges _including_ last one. (histc makes
%extra count for last)

n = histc(x,edges);
n(end-1) = n(end-1) + n(end);
n(end) = [];

end
