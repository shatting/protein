%CLASSSEQLENGTHS main script for viewing and developing the
% used with numsclasseslengthn.m

if (~exist('numscln','var')),
    load RT071127;
    load new_deal;


    sclasses = {[1]};
    %sclasses = {1:27};
    disp('  c_i   c_i+1    t');
    sys10toX(sclasses{1},[3 3 3])

    [ numscln, midx ] = numsclasseslengthn( data, new_deal_data, sclasses );
end

figure;
maxy = max(numscln);
%options
% remove zeros kernel
%zkernelfn = @(x,y) (x+y)/2;
zkernelfn = @(x,y) max(x,y);

% simple kernel
%skernelfn = @(x,y) (x+y)/2;
%skernelfn = @(x,y) max(x,y);
skernelfn = @(x,y) 1/3*min(x,y)+2/3*max(x,y);

% sliding kernel
slkernelfn = @(x) mean(x); %moving average
slkernelfn = @(x) max(x); %moving maximum
%slkernelfn = @(x) sum(sort(x)'.*[1:-1/length(x):1/length(x)])/length(x);
%slkernelfn = @kernel_weighted;

% iterations
nsmooth = 200;
nmatsmooth = 100;
nslsmooth = 1;


% calculate
znum = removezeros(numscln, zkernelfn);

snum = smooth_simple(znum, nsmooth, skernelfn);

slnum = znum;
for i=1:nslsmooth,
    slnum = smooth_sliding(slnum, 10, slkernelfn);
end

combnum = numhelices_filter(numscln);



% plot
if 1, % view combined&final filter
    plot(numscln,'*');
    hold on;
    plot(combnum,'-');
    hold off;
    legend({'original','filtered'});
else % design
    subplot(4,1,1);
    plot(numscln);
    set(gca,'YLim',[0 maxy]);
    title('original');
    
    subplot(4,1,2);
    plot(znum);
    set(gca,'YLim',[0 maxy]);
    title(sprintf('removezeros'));
    
    subplot(4,1,3);
    plot(snum(:,end));
    set(gca,'YLim',[0 maxy]);
    title(sprintf('%i-smoothed',nsmooth));
    
    subplot(4,1,4);
    cfnum = znum;
    for i=1:nmatsmooth,
        cfnum = smooth(cfnum);
    end
    plot(cfnum)
    set(gca,'YLim',[0 maxy]);
    title(sprintf('%i-curve fitting toolbox-smoothed',nmatsmooth));
    
    plot(slnum)
    set(gca,'YLim',[0 maxy]);
    title(sprintf('%i-sliding-smoothed',nslsmooth));
end
