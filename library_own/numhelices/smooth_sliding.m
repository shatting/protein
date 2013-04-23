function [ y ] = smooth_sliding( x, window, kernelfn )
%SMOOTH_SLIDING filter: smooth by sliding windows
%   Detailed explanation goes here
n=length(x);
y = x*0;

if (window == 0),
    for i=1:n,
        y(i) = kernelfn(x(i));
    end
else   
    y(1:window) = kernelfn(x(1:2*window+1));
    y(end-window:end) = kernelfn(x(end-2*window-1:end));

    for i=window+1:n-window-1,
        idx = i-window:i+window;       
        y(i) = kernelfn(x(idx));
    end
end
    
end
