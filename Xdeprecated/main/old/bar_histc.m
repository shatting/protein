function [ h ] = bar_histc( x, y, width, xlims )
%BAR_HISTC draws a stacked bar plot from histc data
% problem is that bar(x,y,'stacked','histc') doesnt work
% fix: shift the domain, because histc uses bin edges whereas usual hist gives
% centers (see histc_test.m)
% assumes fixed width bins

w = x(2) - x(1);

h = bar(x+w/2,y,width,'stacked');

if nargin > 3 && ~isempty(xlims),
    xlim(xlims);
else
    xlim([x(1), x(end)]);
end

end
