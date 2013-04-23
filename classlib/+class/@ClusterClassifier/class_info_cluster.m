function [ class_info ] = class_info_cluster( class_info )
%CLASS_INFO_CLUSTER Summary of this function goes here
%   Detailed explanation goes here
    
    [cluster_info.merge, cluster_info.cost, cluster_info.D0, cluster_info.D] = suff_clusterh(class_info.suff, class_info.options, 1);
    
    c = sum(class_info.suff.count,2);
    cluster_info.freq = squeeze(c(1,:));
    
    cluster_info.readme =  [
    '      merge: see suff_clusterh.m           '
    '       cost: see suff_clusterh.m           '
    '         D0: see suff_clusterh.m           '
    '          D: see suff_clusterh.m           '
    '       freq: class frequencies             '
    '     readme: this readme                   '];    

    class_info.cluster_info = cluster_info;
    
    
    % [ class_info ] = class_info_applyclustering( class_info, relcostbds,
    % numclbds, maxdiff, minprev, minstep )
    class_info = class_info_applyclustering( class_info );
    %class_info_cluster_plot(class_info);
    
end
