function [ class_info ] = class_info_applyclustering( class_info )
%CLASS_INFO_APPLYCLUSTERING
% [ class_info ] = class_info_applyclustering( class_info, relcostbds,
%                                               numclbds, maxdiff, minprev )

if (~isfield(class_info,'cluster_info')),
    disp('no clustering information found, clustering..');
    class_info = class_info_cluster(class_info); 
end

options = class_info.options;

ncl = max(class_info.cl(:,end));
numclbds(1) = ceil(options.clust_minclfact*ncl);
numclbds(2) = ceil(options.clust_maxclfact*ncl);
relcostbds = options.clust_relcostbds;
maxdiff = options.clust_maxdiff;
minprev = options.clust_minprev;
minstep = options.clust_minstep;

% if(isfield(class_info,'cluster_result')),
%    disp('clustering already applied, exiting.');
%    return
% end

cost = class_info.cluster_info.cost;
cost = cost/max(cost);

diff = diffq(cost); % get difference quotients
diff = [-inf diff -inf]; % make same length as cost
nsplits = length(cost);
idx = 1:nsplits;

% xy bounds
numbds = nsplits - numclbds + 1;
numbds = numbds(2:-1:1);
xybds=cost>=relcostbds(1) & cost <= relcostbds(2) & idx>=numbds(1) & idx<=numbds(2);
diff(~xybds)=-inf;

% maxdiff bounds
diffbds = diff<=maxdiff;

steps = cost(2:end)-cost(1:end-1);
steps = [steps 0];
diffbds(steps<minstep) = 0;

% prev bounds
prev = find(cost(2:end)-cost(1:end-1) < minprev);
diffbds(prev+1) = 0;

% apply diffbds
diff(~diffbds) = -1;

[m,cut] = max(diff);

% save some info before possible error
class_info.cluster_info.numbds = numbds;
class_info.cluster_info.numclbds = numclbds;
class_info.cluster_info.diff = diff;
class_info.cluster_info.diffbds = diffbds;

if m<=0,
    warning('clustering not applied, parameters too rigid, please retry later');
else

    merge = class_info.cluster_info.merge(:,1:cut);


    % do merges
    % - into new struct, only cl <-------- doing this
    % OR
    % - merge in suff
    % - get new pot, freqreg, potmin, (cl by simple merges)
    % - update remaining merge to new classes
    % - debug: compare classification by pot with simple merges of cl
    %              pot: [1x1 struct]
    %          freqreg: [19x1 double]
    %           potmin: [3.3317 2.7852 3.4194 3.1600 2.7933 3.5096 2.8326 3.3539 2.6879 1.5301 2.6962 2.5358 2.6098 2.5177 2.7613 2.4757 2.5448 2.5752 2.5829]
    %               cl: [221491x9 uint8]
    %           confus: [3x19 double]
    %             freq: [9x19 double]
    %             suff: [1x1 struct]
    %             type: [400 400 400 400 400 400]
    %          options: [1x1 struct]
    %           readme: [12x43 char]
    %          featfun: @feat_pairs
    %             time: 18.9996
    %     cluster_info: [1x1 struct]

    cl = applymerge(class_info.cl(:,end), merge);

    % save results
    cluster_result.cl = cl;
    cluster_result.ncl = max(cl);
    cluster_result.cut = cut;

    class_info.cluster_result = cluster_result;

    cluster_result
end

end