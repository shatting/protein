function tplot_08( ts, cl, tit, xlab, ylab )
%TPLOT_08 plots 

n = size(ts,1);

ts = (double(ts) + rand(n,1))/100 * pi;

p = [cos(ts) sin(ts)];

rr = 0.2;
r = rand(n,1)*rr + 1-rr;

%plot(r.*p(:,1),r.*p(:,2),'.')
%nf=gcf;
%figure(99)
%hist(ts,200);
%input('disc>');
%figure(nf)
scatterplot(p.*[r r],cl,'','',n,1);
title(tit);
xlabel(xlab);
ylabel(ylab);
hold on;

r = 0.55;
for i=-90:10:100,
    ang = (i)/100 * pi;

    text(r*cos(ang)+0.45, r*sin(ang)+0.5, int2str(i));
end

%axes off;

hold off;
set(gca,'xtick',[])
set(gca,'ytick',[])
axis square;

end
