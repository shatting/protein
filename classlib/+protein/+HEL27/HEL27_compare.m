function [ ] = new_deal_compare( nd1, nd2 )
%NEW_DEAL_COMPARE 

if (nargin==0 || ~exist([nd1,'.mat'],'file') || ~exist([nd2,'.mat'],'file')),
    dir('new_deal*.mat');
    return
end

load([nd1,'.mat']);
ndd1 = new_deal_data;
dprintf('--------------------------\n%s.mat',nd1)
dprintf('pre-run comment:\n\t\t%s',ndd1.comments);
dprintf('post-run comment:\n\t\t%s',ndd1.pr_comments);
dprintf('clustered: %i',isfield(ndd1.hefg_info,'cluster_result'));
load([nd2,'.mat']);
ndd2 = new_deal_data;
dprintf('--------------------------\n%s.mat',nd2)
dprintf('pre-run comment:\n\t\t%s',ndd2.comments);
dprintf('post-run comment:\n\t\t%s',ndd2.pr_comments);
dprintf('clustered: %i',isfield(ndd2.hefg_info,'cluster_result'));

dprintf('\n\n---> from now on: rows/first class set=%s, columns/second class set=%s <----\n\n',nd1,nd2);

disp('++++++++ NON-CLUSTERED classes +++++++++++');
cl1 = ndd1.hefg_info.cl(:,end);
cl2 = ndd2.hefg_info.cl(:,end);

new_deal_quality(cl1,cl2,1,80,1);

disp('++++++++++++++++++++++++++++++++++++++++++');
if (isfield(ndd1.hefg_info,'cluster_result') && isfield(ndd2.hefg_info,'cluster_result')),
    disp('++++++++ CLUSTERED classes +++++++++++');
    cl1 = ndd1.hefg_info.cluster_result.cl;
    cl2 = ndd2.hefg_info.cluster_result.cl;

    new_deal_quality(cl1,cl2,1,80,1);
    disp('++++++++++++++++++++++++++++++++++++++++++');
end

end
