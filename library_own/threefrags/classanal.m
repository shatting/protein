close all;
potcary = [potc{:}];
means = [potcary(:).mean];

[U,S,V] = svd(means);

dprintf('singular values of means:');
diag(S)'

clend = cl(:,end);
clnonempty(nonempty) = clend; % get raw classes
ncl = max(clend);

nrbins = length(rawdata.rbins)-1

centeraas = cell(ncl,1);
suffgrps = centeraas;
plotx = 4;
ploty = 4;

for c=1:ncl,
    clidx = find(clend==c);
    aarbin = sys10toX(clidx,[20 20 nrbins]);
    
    suffgrps{c} = [suffdata.suffs{nonempty(clidx)}];
    
%     % put group mean r into 4th column
%     for i=1:length(suffgrps{c}),
%         s = suffgrps{c}(i);
%         aarbin(i,4) = s.moment(end,11)/s.moment(end,end);
%     end
    
    centeraas{c} = aarbin;
    fig = floor((c-1)/(plotx*ploty));
    figure(100+fig)
    set(gcf,'name',sprintf('histograms page classes %i to %i',fig*plotx*ploty+1,(fig+1)*plotx*ploty));
    subplot(ploty,plotx,mod(c-1,plotx*ploty)+1);
    hist(aarbin,20);
    title(sprintf('class %i: %i groups, %i points',c,length(clidx),potc{c}.freq));
    %legend({'i-aa','j-aa','rbin'})
end

figure;
set(gcf,'name','aa combination distribution per class');
x = ceil(sqrt(ncl));
for c=1:ncl,
    clidx = find(clend==c);
    aarbin = sys10toX(clidx,[20 20 length(rawdata.rbins)-1]);
    
    subplot(x,x,c);
    spy(sparse(aarbin(:,1),aarbin(:,2),aarbin(:,3),20,20))
    xlabel '';
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
       
end;

figure;
set(gcf,'name','bin distributions per class');
x = ceil(sqrt(ncl));
for c=1:ncl,
    clidx = find(clend==c);
    aarbin = sys10toX(clidx,[20 20 nrbins]);
    
    subplot(x,x,c);
    spy(sparse(aarbin(:,3),aarbin(:,1),aarbin(:,3),nrbins,20))
    xlabel '';
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
       
end;

