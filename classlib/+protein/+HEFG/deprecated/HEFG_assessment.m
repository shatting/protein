close all;

hefg_getfeats;

hefg = double(new_deal_data.hefg_info.cl(:,1));
ts = feats(:,3:4);
cs = feats(:,1:2);
seq = double(feats(:,5:8));


if (~isfield(new_deal_data.hefg_info,'cluster_result') || ask('clustering info found. assess non-clustered classes?')),
    cl = double(new_deal_data.hefg_info.cl(:,end));
    clso = sortcl(cl);
    classification_assessment([hefg clso], cs, ts, seq, 10)
    new_deal_quality(hefg, clso, 1,90);
end

if (isfield(new_deal_data.hefg_info,'cluster_result') && ask('clustering info found. assess cluster result classes?')),
    class_info_cluster_plot(new_deal_data.hefg_info);
    
    cl = double(new_deal_data.hefg_info.cluster_result.cl);
    clso = sortcl(cl);
    classification_assessment([hefg clso], cs, ts, seq, 10);
    new_deal_quality(hefg, cl, 1,90);
end
