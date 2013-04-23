function HEL27_info( matfilename )
%NEW_DEAL_INFO Summary of this function goes here
%   Detailed explanation goes here

global classdatadir;
load([classdatadir,filesep,matfilename]);

whos

c_class_info = new_deal_data.c;
cp_class_info = new_deal_data.cp;
final_class_info = new_deal_data.final;
comments = new_deal_data.comments;
pr_comments = new_deal_data.pr_comments;

disp('_______________________ c_i _______________________________');
class_info(c_class_info)
if (isfield(c_class_info,'cluster_info'))
    %c_class_info.cluster_info
    cluster_info_plot(c_class_info.cluster_info, c_class_info.cluster_result);
    title 'c_i clusters';
end

disp('_______________________ c_i+1 _______________________________');
class_info(cp_class_info)
if (isfield(cp_class_info,'cluster_info'))
    %cp_class_info.cluster_info
    cluster_info_plot(cp_class_info.cluster_info, cp_class_info.cluster_result);
    title 'c_{i+1} clusters';
end

disp('_______________________ final _______________________________');
class_info(final_class_info)
if (isfield(final_class_info,'cluster_info'))
    %final_class_info.cluster_info
    cluster_info_plot(final_class_info.cluster_info, final_class_info.cluster_result);
    title 'final clusters';
end

disp('______________________________________________________');
pre_run_comments = comments

post_run_comments = pr_comments

end
