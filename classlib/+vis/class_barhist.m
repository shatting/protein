function [] = class_barhist( data, cl, nbins, smear, nbiggest, stretchy, barwidth, inlabel, stretchx)
%class_barhist( data, cl, nbins, integerdata, nbiggest, stretchy, barwidth, inlabel )
% produces a 'histogram' of data(:,1) with classe cl
%   nbins       .. either number of bins or, for integer data n, [n_min:1:n_max]
%   smear       .. if yes, random quantity \in [-1,1] added to data
%   nbiggest    .. draw only nbiggest biggest classes 
%   stretchy    .. stretch bins vertically, makes histogram 'relative'
%   barwidth    .. width of bars in (0,1], default 1 if omitted
%   inlabel     .. no legend if [], legend with inlabel labels if cell array, or class
%                  numbers if omitted
%   stretchx    .. if 1, data is transformed s.t. all bins have the same
%                  number of samples (makes sense only with stretchy==1 and no integer data)

data = double(data(:,1));
mind = min(data);
maxd = max(data);
cl = double(cl);
n = size(data,1);
ncl = max(cl);

if (nargin<6), stretchy = 0;end
if (nargin<9), stretchx = 0;end

if nargin >= 4 && smear, 
    data = min(max(data + rand(size(data,1),1)*2 - 1,mind),maxd);
end

if nargin >= 5 && nbiggest>0, 
    stopidx = max(ncl - nbiggest + 1,1);
else
    stopidx = 1;
end

if (nargin<7)
    barwidth = 1;
end

map = class_colormap(ncl);
freq = freqs(cl)/n;

if length(nbins)==1,
    domain = (maxd-mind);
    binwidth = domain/nbins;
    edges = mind:binwidth:maxd;    
    
    
    if (stretchx),
         ed = histc_fix(data,edges);
                 
         [s,i] = sort(data);
         perbin = floor(n/nbins);
         lastbin = n - perbin*nbins;
         
         for i=1:nbins,
            ix = (i-1)*perbin;
            xedges(i) = s(ix+1);
         end
         
         xedges(nbins+1) = maxd;
    end
    xlims = [];
else
    edges = [0 nbins] + 0.5;  % integers
    nbins = max(nbins);
    xlims = [0.5 nbins+0.5];
end
xc = zeros(nbins,nbiggest);

%sort descending by size for drawing largest first
[x,asortcl] = sort(freqs(cl));
labl={};

%get total counts
for i=ncl:-1:stopidx,
    iplot = ncl-i+1;
    c = asortcl(i);
    if nargin < 8,
        labl{iplot} = int2str(c);
    elseif ~isempty(inlabel)
        labl{iplot} = inlabel{c};
    %else
        % labl remains empty
    end
    
    cdata = double(data(cl==c));
    
    if (~stretchx)
        cc = histc_fix(cdata,edges);
    else
        cc = histc_fix(cdata,xedges);
    end
    
    xc(:,iplot) = cc;
end

if (stretchy)
    xc = xc./repmat(sum(xc,2),1,size(xc,2)); % make rows sum to 1
    ylim = [0 1];
%    set(gca,'XAxis','off');
else
    ylim = [0 max(sum(xc,2))];
end

maximize_fig;

bh=bar_histc(edges(1:end-1),xc,barwidth,xlims);  % use x range from nlc-th centers (should all be the same anyway..), width 1


colormap(map(asortcl(ncl:-1:stopidx),:));
if ~isempty(labl)
    lh=legend(labl);
    y = get(lh,'outerposition');
    y(1) = 0.91;
    set(lh,'outerposition',y);
end

if nargin>5 && stretchy, set(gca,'YTick',[]); end

%set(gca,'XLim',[mind maxd]);

set(gca,'YLim',ylim);

set(gca,'box','off','tickdir','out','xminortick','on','fontsize',12);


end
