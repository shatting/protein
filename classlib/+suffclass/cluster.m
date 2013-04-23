%SUFF_CLUSTER Cluster analysis, cut finding and applying.

function [ cl ] = suff_cluster( suff, cl, options )
   
    [merge, cost, D0, D] = suff_clusterh(suff, options, 1);
    
    c = sum(suff.count,2);
    freq = squeeze(c(1,:));
    
    cut = cluster_tryfindcut(cl,cost,options);
     % truncate merge at cut
    merge = merge(:,1:cut);

    % do merges
    % - into new struct, only cl <-------- doing this
    % OR
    % - merge in suff
    % - get new pot, freqreg, potmin, (cl by simple merges)
    % - update remaining merge to new classes
    % - debug: compare classification by pot with simple merges of cl

    cl = applymerge(cl, merge);

    % TODO: visualize results
        
end
