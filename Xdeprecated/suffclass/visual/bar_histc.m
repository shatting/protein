%BAR_HISTC Draws stacked bar plot from histc data.
% problem is that bar(x,y,'stacked','histc') is not implemented in MATLAB.
% fix: shift the domain, because histc uses bin edges whereas usual hist gives
% centers
% assumes fixed width bins
%
% See also: bar_histc_test.


function [ h ] = bar_histc( x, y, width, xlims )

w = x(2) - x(1);

h = bar(x+w/2,y,width,'stacked');

if nargin > 3 && ~isempty(xlims),
    xlim(xlims);
else
    xlim([x(1), x(end)]);
end

end
