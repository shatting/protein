function [ s ] = kldiv_best_classifier(verbose)
%KLDIV_BEST_CLASSIFIER find best hefg classifier in directory

x = dir('new_deal_hefg*.mat');

scores = zeros(length(x),2) + inf; %[normal clustered]

sizes = zeros(length(x),2);

for i=1:length(x),
    load(x(i).name);
    
    isclust = isfield(new_deal_data.hefg_info,'cluster_result');
    
    dprintf(' ++++++++++++++++ %s ++++++++++++++++++++ ',x(i).name);
    dprintf('pre-run comment:\n\t\t%s',new_deal_data.comments);
    dprintf('post-run comment:\n\t\t%s',new_deal_data.pr_comments);
    dprintf('clustered: %i',isclust);

    cl_hefg = double(new_deal_data.hefg_info.cl(:,1));
    cl_nocl = double(new_deal_data.hefg_info.cl(:,end));    
    sizes(i,1) = max(cl_nocl);
    
    dprintf('unclustered ------------------------------------------');
    dprintf('number of classes: %i',max(cl_nocl));
    [q] = kldiv2qualitymeasure(cl_hefg, cl_nocl,'sym');
    dprintf('final measure: %f',q);
    scores(i,1) = q;
    
    if (isclust)
        dprintf('clustered   ------------------------------------------');
        cl_cl = double(new_deal_data.hefg_info.cluster_result.cl);
        sizes(i,2) = max(cl_cl);
        dprintf('number of classes: %i',max(cl_cl));
        [q] = kldiv2qualitymeasure(cl_hefg, cl_cl,'sym');
        scores(i,2) = q;
        dprintf('final measure: %f',q);
    end
    
end

dprintf('+++++++++++++++++++++++++++++++++++++++++++++++\nSUMMARY.');
dprintf('normal          clustered       name  ');        

%scores = scores./sizes;

for i=1:size(scores,1),
    dprintf('%4.4f (%i)  \t|  %4.4f (%i)  \t| %s',scores(i,1),sizes(i,1),scores(i,2),sizes(i,2),x(i).name);
end

[m,i] = min(scores);
[k,l] = min(m);

mindeal = i(l);
clustered = l;
dprintf('+++++++++++++++++++++++++++++++++++++++++++++++\n...and the WINNER is...');
dprintf('%4.4f (%i) | clustered = %i | %s',scores(mindeal,clustered),sizes(mindeal,clustered),clustered-1,x(mindeal).name);


end
