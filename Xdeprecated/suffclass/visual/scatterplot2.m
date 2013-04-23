function [] = scatterplot2( data, data2, idxknown, mu, C, suff)
%SCATTERPLOT2 used in suff_recon_test

[n,d] = size(data);

if nargin > 5,
    pot = suff_num2pot(suff,0);
end
figure;

for i = 1:d-1,
    for j = i+1:d,
        subplot(d-1,d-1,(d-1)*(i-1)+j-1);

        plot(data(:,i),data(:,j),'k.');
        hold on;
        plot(data2(:,i),data2(:,j),'ro');
        
        %ellipse(mu([i,j]),1,schurComp(C,[i,j]),'k-');
        ellipse(mu([i,j]),1,C([i,j],[i,j]),'k-');
        if nargin > 5, ellipse(pot.mean([i,j]),2,pot.cov([i,j],[i,j]),'b.-'); end
        %ellipse(mu,1,C(1:2,1:2),'k:');
        %ellipse(pot.mean,1,pot.cov(1:2,1:2),'r:');
        set(gca,'XTick',[]);
        set(gca,'XTickLabel',[]);
        set(gca,'YTick',[]);
        set(gca,'YTickLabel',[]);        

        if ismember(i,idxknown) si = '*'; else si = ''; end
        if ismember(j,idxknown) sj = '*'; else sj = ''; end        
        
        xlabel(sprintf('%i%s',i,si));
        ylabel(sprintf('%i%s',j,sj));        
        %title(sprintf('%i%s/%i%s',i,si,j,sj));
        
        if ismember(i,idxknown) || ismember(j,idxknown)
            set(gca,'color','cyan');
        end
        if ismember(i,idxknown) && ismember(j,idxknown)
            set(gca,'color','green');
        end
        axis equal;
        hold off;
    end
end
legend1 = legend({'original data','reconstruced data','true covariance','estimated covariance'});
set(legend1,'Orientation','horizontal','Location','NorthWestOutside');

end