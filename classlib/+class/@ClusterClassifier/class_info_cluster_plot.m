function [ g2c ] = class_info_cluster_plot( class_info )
%CLUSTER_INFO_PLOT 
% cluster_info_plot( cluster_info, cluster_result, target, cl )
% cluster_info_plot( class_info )


    %disp('merge:');
    %tightmat(cluster_info.merge);
    cluster_info = class_info.cluster_info;
    options = class_info.options;    
    
    figure;
    % --------------
    subplot(3,2,[1,3,5]);
    g2c = dendro(cluster_info.freq, cluster_info.merge, cluster_info.cost, cluster_info.D0);

    % --------------
    subplot(3,2,[2,4]);
    hold on;
    xr = [0,length(cluster_info.cost)+1];
    yr = [0,101];
    
    cost = cluster_info.cost*100/max(cluster_info.cost);


    idxok = find(cluster_info.diff>0);
    idxnok = find(cluster_info.diff<=0);
    % plot ok cuts
    plot(idxok,cost(idxok),'+g');    

    % plot not ok cuts
    plot(idxnok,cost(idxnok),'+r');    

    if isfield(class_info,'cluster_result'),
        % plot chosen cut
        plot(class_info.cluster_result.cut,cost(class_info.cluster_result.cut),'ob');
    end
        
    % plot quality
    nmerge = length(cost);
    q = zeros(1,nmerge);
    clm = class_info.cl(:,end);
    target = class_info.cl(:,1);
    mg = cluster_info.merge;
    % successively merge
    for i=1:nmerge-1,
        %dprintf('merge %i',i);
        [clm, rep] = applymerge(clm,mg(:,1));
        % rename classes in merge
        mg = mg(:,2:end);
        mg = rep(mg);
        qi = kldivqualitymeasure(target, clm, 0);
        q(i) = qi;
    end
    plot(q*100,'.b');

    % plot xy bounds
    plot(xr,[options.clust_relcostbds(1) options.clust_relcostbds(1)]*100-0.5,'-.r');
    plot(xr,[options.clust_relcostbds(2) options.clust_relcostbds(2)]*100+0.5,'-.r');
    plot([cluster_info.numbds(1) cluster_info.numbds(1)]-0.5,yr,'-.r');
    plot([cluster_info.numbds(2) cluster_info.numbds(2)]+0.5,yr,'-.r');
    hold off;
    xlim(xr);
    ylim(yr);

    % --------------
    subplot(3,2,[6]);
    hold on;
    diff = diffq(cluster_info.cost);
    diff = [0 diff 0];
    diff(diff>options.clust_maxdiff) = options.clust_maxdiff;
    
    idxok = find(cluster_info.diffbds);
    idxnok = find(~cluster_info.diffbds);
        
    plot(idxok,diff(idxok),'+g');

    plot(idxnok,diff(idxnok),'+r');
               
    if isfield(class_info,'cluster_result'),
         plot(class_info.cluster_result.cut,diff(class_info.cluster_result.cut),'ob');
    end
    
    hold off;
    xlim([0,length(cluster_info.cost)+1]);
    ylim([0 options.clust_maxdiff]);   

end
