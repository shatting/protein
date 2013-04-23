function [] = class_hist( data, cl, suby, subx, subidx, nbins, xlim, ylim, integerdata, nbiggest )
%class_hist( data, cl, suby, subx, subidx, nbins, xlim, ylim, integerdata, nbiggest )
%   subx, suby ... width, height of subplot
%   subidx     ... subplot indexes [idx1,..,idxn]
%   integerdata .. is data integerdata. if yes, random quantity \in [-1,1] added

data = double(data);
mind = min(data);
maxd = max(data);
cl = double(cl);

n = size(data,1);
ncl = max(cl);

if nargin >= 9 && integerdata, 
    data = min(max(data + rand(size(data,1),1)*2 - 1,mind),maxd);
end

if nargin >= 10 && nbiggest>0, 
    stopidx = ncl - nbiggest + 1;
else
    stopidx = 1;
end


map = class_colormap(ncl);
freq = freqs(cl)/n;
nbinsc = round(freq*nbins);

if (nargin>6 && ~isempty(xlim))
    
else
    xlim = [mind,maxd];
end

step = (xlim(2) - xlim(1))/(nbins+1);
edges = xlim(1):step:xlim(2);

%sort descending by size for drawing largest first
[x,asortcl] = sort(freqs(cl));
curmax = -1;

for i=ncl:-1:stopidx,
    c = asortcl(i);
    labl{ncl-i+1} = int2str(c); % label    
    cdata = double(data(cl==c)); % include minimum and maximum so that the range is the same for all classes implying same bins over all classes
    
    subplot(suby,subx,subidx);
        
    clct = histc(cdata,edges);
    
    bar(edges,clct,'histc')

    curmax = max([clct;curmax]);
    hold on
    h = findall(gca,'Type','patch','-not','Tag','colored');
    set(h,'Tag','colored');
    set(h,'FaceColor',map(c,:));
    set(h,'EdgeColor',map(c,:));
end
hold off

legend(labl);

if ~isempty(xlim),
    set(gca,'XLim',xlim);
end
if nargin >= 8 && ~isempty(ylim),    
    set(gca,'YLim',ylim);
else
    set(gca,'YLim',[0 ceil(curmax*1.1)]);
end

set(gca,'box','off','tickdir','out','xminortick','on','fontsize',12);

end
