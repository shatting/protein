function [ outp ] = numhelices_filter( inp )
% combined&final filter
% 1. remove zeros with dynamic window and maximum kernel
% 2. do a fixed window moving maximum
% 3. make it look nice with a fixed window moving average
% 4 finally, ceil() all values
combfilt = @(x) smooth_sliding(smooth_sliding(smooth_sliding(removezeros(x, @(x,y) max(x,y)) ,10, @(x) max(x)), 15, @(x)mean(x)),0,@(x)ceil(x));
outp = combfilt(inp);

end
