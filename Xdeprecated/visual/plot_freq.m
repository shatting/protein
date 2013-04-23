function [] = plot_freq( freq )
%FREQ_PLOT 

ncl = size(freq,2);
numit = size(freq,1);
colmap = class_colormap(ncl);

figure;
hold on;

for c = 1:ncl,
   plot(freq(:,c),'Color',colmap(c,:));
end

hold off;

xlabel('iteration');
ylabel('frequency');
legend ('show');
set(gca,'XTick',1:numit);