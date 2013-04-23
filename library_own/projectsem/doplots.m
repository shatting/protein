
if ~exist('histogram'),
    threefragments;
end

for aa3=1:20,%size(histang,3),   
    disp(sprintf('printing sheet %i',aa3));
    figure(aa3);
    clf;    
    hold on;    
    xlabel('aa1');
    ylabel('aa2');
    title(sprintf('(aa2,aa1,%i)', aa3));
    for aa1=1:20,
        disp(sprintf('row %i',aa1));
        for aa2= 1:20,
               
               histo(histogram(aa1,aa2,aa3,:),[aa1,aa2], sprintf('%i,%i,%i [%i]',aa1,aa2,aa3,sum(histogram(aa1,aa2,aa3,:))));

        end
    end
        
    colormap('hot');
    axis off;
    hold off;
    saveas(gcf,sprintf('seite_%i.jpg',aa3));
end

