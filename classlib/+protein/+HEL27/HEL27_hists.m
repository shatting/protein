close all

[HEL27, HEL27_ccpt] = HEL27_anal08(geom);

%HEL27_ccpt = sys10toX(HEL27,[3 3 3]); % c_i, c_i+1 and t classes

figure('Name','c_i and c_{i+1} - dependent histograms for t');

class_hist(geom(:,1), HEL27_ccpt(:,1), 4, 3, [1 2 3], 160, [-60 100], [0 10000]);

xlabel(sprintf('c, HEL==rgb (n=[%i,%i,%i])',sum(HEL27_ccpt(:,1)==1),sum(HEL27_ccpt(:,1)==2),sum(HEL27_ccpt(:,1)==3)));


for ci=1:3,
    for cip=1:3,
        indci = HEL27_ccpt(:,1) == ci;
        indcip = HEL27_ccpt(:,2) == cip;
        ind = indci & indcip;
        geoms = geom(ind,:);
        n = size(geoms,1);
        if 1 % training
          HEL12=[ci,cip]
          figure(ci+3*cip-2);
          subplot(2,1,1);
          class_hist(geom(ind,3),HEL27_ccpt(ind,3),2,1,1,2000,[-100,100]);
          title(sprintf('c_i = %i, c_{i+1}=%i',ci,cip));
          subplot(2,1,2);
          tplot_08(geoms(:,3),HEL27_ccpt(ind,3),'',sprintf('n=%i',n),sprintf('c_i = %i, c_{i+1} = %i',ci,cip));
          input('next>');
        end;
        figure(1)
        subplot(4,3,ci+3*cip);
        tplot_08(geoms(:,3),HEL27_ccpt(ind,3),'',sprintf('n=%i',n),sprintf('c_i = %i, c_{i+1} = %i',ci,cip));
    end
end

