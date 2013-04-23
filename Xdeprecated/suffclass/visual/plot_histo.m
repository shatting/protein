%PLOT_HISTO Plot several histograms in one figure.
%   data   histogram bins
%   pos    position in grid, eg [3,8]
%
% See also: plot_histo_test.

function [  ] = plot_histo( data, pos, labeltext )

XWIDTH = 100;
YHEIGTH = 80;

nbins = length(data);
binwidth = XWIDTH/nbins;

count = sum(data);
maxcount = max(data);

xoff = pos(1)*XWIDTH*1.2;
yoff = -pos(2)*YHEIGTH*1.3;

if (maxcount>0),
    for b=1:nbins,
        x = [xoff+binwidth*b, xoff+binwidth*b, xoff+binwidth*(b+1), xoff+binwidth*(b+1)];
        y = [yoff, yoff+YHEIGTH*data(b)/maxcount, yoff+YHEIGTH*data(b)/maxcount, yoff];
        fill(x,y,-count);
        
        %plot([xoff+binwidth*b, xoff+binwidth*b],[yoff, yoff+YHEIGTH*data(b)/maxcount]);
        %plot([xoff+binwidth*b, xoff+binwidth*(b+1)],[yoff, yoff]);
        %plot([xoff+binwidth*b, xoff+binwidth*(b+1)],[yoff+YHEIGTH*data(b)/maxcount, yoff+YHEIGTH*data(b)/maxcount]);
        %plot([xoff+binwidth*(b+1), xoff+binwidth*(b+1)],[yoff, yoff+YHEIGTH*data(b)/maxcount]);
    end;    
end

%text(xoff+XWIDTH/2,yoff+YHEIGTH,labeltext);
%set(h,'fontsize',6);
