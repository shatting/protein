function [] = plot_t( p, cl, numcols, titl )
%TPLOT plot t vs t' using numcols colors
% plot_t( p, cl, numcols, titl )
% some errors involving class sizes (27.6.)

figure;
maximize_fig;

if nargin < 3,
    numcols = 2;
end
if any(cl == 0),
    cl(cl == 0) = max(cl) +1;
end

if (isempty(cl))
   cl = ones(size(p,1),1); 
end

classes = unique(cl);
numclasses = length(classes);
numrows = ceil(numclasses / numcols) + 2;

freq = freqs(cl);
[scl, sidx] = sort(freq);
invp(numclasses:-1:1) = sidx;
invp(invp) = 1:numclasses;

lefts = floor(numcols/2);
lindex = [1:lefts,numcols+1:numcols+lefts];
rindex = setdiff([1:2*numcols],lindex);

subplot(numrows,numcols,lindex);
scatterplot(p,invp(cl),'','',100000, 0.1, 0);
axis square;
xlims = [min(p(:,1)) max(p(:,1))];
ylims = [min(p(:,2)) max(p(:,2))];
title(titl);

subplot(numrows,numcols,rindex);
%bar(freq);
class_barhist(cl,cl,max(cl),0,max(cl),0);


for c=1:numclasses,   
    subplot(numrows, numcols, c + 2*numcols);
    clsnum = sidx(numclasses - c + 1);
    scatterplot(p(cl == clsnum,:),clsnum,'','',10000, 0.1, 0, 0, nclasses);
    title(sprintf('class %i, n=%i (%2.0f%%)',clsnum, freq(clsnum),100*freq(clsnum)/length(cl)));
    xlim(xlims);
    ylim(ylims);
    axis square;
end;