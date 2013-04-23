option.newdata = 1;
option.filter = 0;

if option.newdata,
    addpath '../suffclass';
    addpath ../utils;

    tic
    close all;

    % setup
    classsigma = [1.5,1.5,1.5];
    classmu = [[-1.8,-1.8];[0,2];[1,-3]];
    oneclassparam = 0.75;
    twoclassparam = 0.39;
    threeclassparam = 0.28;
    numgroups = 700;
    numseedperclass = 48; % so we hopefully get around 12 in the initial classification
    numseedclasses = 7;


    numclasses = length(classsigma);
    numpoints = fix(max((randn(numgroups + numseedperclass*numseedclasses,1) + 3) * 200, 20)); % numper of points in group

    figure;
    hist(numpoints);
    title 'group density';

    groupdist = rand(numgroups,numclasses);                   % relative distribution of classes in group
    tot = sum(groupdist,2);
    groupdist = groupdist ./ repmat(tot,1,numclasses);

    seeddensities(1,:,:,:) = [0.9,0.05,0.05];
    seeddensities(2,:,:,:) = [0.05,0.9,0.05];
    seeddensities(3,:,:,:) = [0.05,0.05,0.9];
    seeddensities(4,:,:,:) = [0.45,0.45,0.1];
    seeddensities(5,:,:,:) = [0.1,0.45,0.45];
    seeddensities(6,:,:,:) = [0.45,0.1,0.45];
    seeddensities(7,:,:,:) = [0.33,0.33,0.34];

    % add seeds to groups
    groupdist = [groupdist;repmat(seeddensities,numseedperclass,1)];

    groupdistabs = fix(groupdist .* repmat(numpoints,1,numclasses));

    totalnumpoints = sum(sum(groupdistabs));

    % generate covariances for the classes
    for i=1:numclasses,    
        t = zeros(2,2);
        while cond(t) > 8;
            t = (rand(2,2)-0.5);   
            t = t*t'*classsigma(i);
        end
        classcovs{i} = t;
    end  

    % generate group seed densities
    p = zeros(totalnumpoints,2);

    %generate group points
    figure;
    clf;
    hold on;
    colmap = colormap(hsv(numclasses));
    groupmap = colormap(hsv(numgroups));

    groupdata = cell(numgroups,1);

    index = 1;
    for i=1:numgroups,

        pg = [];
        for j = 1:numclasses,
            pc = genpoints(classmu(j,:), classcovs{j}, groupdistabs(i,j));

            % plot points
            if (rand > 0.1)
                plot(pc(:,1),pc(:,2),'.','Color',colmap(j,:));       
            end
            pg = [pg; pc];
        end

        idx = index:index+size(pg,1)-1;
        p(idx,:) = pg;
        groupdata{i} = pg;

        index = index + size(pg,1);
        %plot group
        %plot(pg(:,1), pg(:,2), 'o', 'Color', groupmap(i,:));

    end

    xlims = [min(p(:,1)) max(p(:,1))];
    ylims = [min(p(:,2)) max(p(:,2))];

    % plot class mu and cov
    for i = 1:numclasses,
       plot(classmu(i,1), classmu(i,2),'x','Color',colmap(i,:));
       ellipse(classmu(i,:),1,classcovs{i},'-',colmap(i,:));
    end

    title 'raw data';
    hold off;


    % create suff stats
    suffs = cell(numgroups,1);

    emptystat = suffstat(1,zeros(2,1));

    index = 1;
    for g=1:numgroups,
    %   figure(100 + j);
    %   hold on;
        data = groupdata{g};
        suffs{g} = suffstat(emptystat, data, size(data,1), ones(size(data,1),1));

    %   pot = suffpot(suffs{index});      
    %   plot(grouppoints(:,1),grouppoints(:,2),'.');        
    %   ellipse(pot.mean,1,pot.cov,'-k');
    %   hold off;
    end

    initclasses = zeros(numgroups,1);
    % set some initial classes
    for i=1:numgroups,
        if (groupdist(i,1) > oneclassparam)
            initclasses(i) = 1;
        end
        if (groupdist(i,2) > oneclassparam)
            initclasses(i) = 2; 
        end
        if (groupdist(i,3) > oneclassparam)
            initclasses(i) = 3; 
        end
        if (groupdist(i,1) > twoclassparam && groupdist(i,2) > twoclassparam)
            initclasses(i) = 4; 
        end
        if (groupdist(i,2) > twoclassparam && groupdist(i,3) > twoclassparam)
            initclasses(i) = 5; 
        end
        if (groupdist(i,1) > twoclassparam && groupdist(i,3) > twoclassparam)
            initclasses(i) = 6; 
        end
        if (groupdist(i,1) > threeclassparam && groupdist(i,2) > threeclassparam && groupdist(i,3) > threeclassparam)
            initclasses(i) = 7; 
        end
    end

    % plot initial assignments
    suffclass_test_plot2(0,suffs,initclasses,0,1);
    title 'seed points + group covariance';

    figure;
    hist(initclasses,7);
    title 'histogram of initial classes';

end % newdata

% filter
if option.filter,
    filterindex = initclasses<4;
else
    filterindex = 1:length(initclasses);
end

filteredsuffs = suffs(filterindex);
filteredinitclasses = initclasses(filterindex);

% classify
[pot,cl,confus, freq] = suffcov(filteredsuffs,filteredinitclasses);

% plot class assignment for raw data
numfinalclasses = max(cl);
classmap = colormap(hsv(numfinalclasses));
for c=1:numfinalclasses,
    idx_final = [1:length(cl)];
    idx_final = idx_final(cl' == c);
    if (length(idx_final) == 0)
        continue
    end
    
    figure;

    d = groupdata(idx_final);
    a = [];
    for i=1:length(d),
        a = [a;d{i}];
    end
        
    subplot(2,3,1);
    hold on;
    plot(a(:,1), a(:,2), '.', 'Color', classmap(c,:));
    xlim(xlims); ylim(ylims);    
    finalnum = sum(cl==c);
    title(sprintf('raw data for final class %i (n=%i, %2.0f%%)',c, finalnum, finalnum/length(initclasses)*100));
    hold off;
    
    idx = [2,3,5,6];
    k = 1;
    for i=idx,
        if k > length(d),
            break
        end
        subplot(2,3,i);
        x = d{k};
        plot(x(:,1), x(:,2), '.', 'Color', classmap(c,:));        
        xlim(xlims); ylim(ylims);
        g = groupdist(idx_final(k),:);
        title(sprintf('example group %i (%2.0f%%, %2.0f%%, %2.0f%%)',k,g(1)*100,g(2)*100,g(3)*100));
        k = k + 1;
    end

    d = groupdata(filteredinitclasses == c);
    a = [];
    for i=1:length(d),
        a = [a;d{i}];
    end
    subplot(2,3,4);
    hold on;
    plot(a(:,1), a(:,2), '.', 'Color', classmap(c,:));
    xlim(xlims); ylim(ylims);
    initnum = sum(initclasses==c);
    
    title(sprintf('raw data for initial class %i (n=%i, %2.0f%%)',c, initnum, initnum/length(initclasses)*100));
    hold off;
end

figure;
hold on;
for c=1:numfinalclasses,  
    idx = cl' == c;
    if (max(idx) == 0)
        continue
    end
    d = groupdata(idx);
    a = [];
    for i=1:length(d),
        a = [a;d{i}];
    end
    plot(a(:,1), a(:,2), '.', 'Color', classmap(c,:));    
end
title(sprintf('raw data for all final classes'));
hold off;

plot_freq(freq(:,1:7))
title('class frequencies');
