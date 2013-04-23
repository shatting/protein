%% HEL27
%% data loading

clc;
close all;
diary off;
if (exist('HEL27.dry','file')), delete('HEL27.dry'); end;

sclassifier = class.HEL27SClassifier();
hel27 = sclassifier.classify(data.FeatureDB.db);
fragaaseq = data.FeatureDB.db.getdata({'aa1','aa2','aa3','aa4'});

comments = sprintf(input('pre-run comments: ','s'));


usediary = 0;
if (usediary)                
    diary HEL27.dry
end

% options
clustercs = 0;
clustert = 0;
mintriples = 300;
numtriples = 100; % override mintriples if > 0
featfun = @feat_pairs;
%featfun = @(x) x;

options = struct(   'plots',0,...
                    'max_iterations',1,...
                    'ask',0,...
                    'ask_timeout',5,...
                    'term_percentage',0.05,...
                    'new_classes',0.015,...
                    'stop_on_nonewclasses',0,...
                    'num_reg_weight',1,...
                    'num_ent_weight',1,...
                    'cat_prior_reg',0.01,...
                    'cat_cond_reg',realmin,...
                    'cat_ent_weight',0.5,...
                    'cat_no_lookup',0,...
                    'final_adjust',0,...
                    'savedata',0);
%% data preparation                
fseq = featfun(fragaaseq);
HEL27_ccpt = sys10toX(hel27,[3 3 3]);

%% c classification
tic;
% find prediction for c_i     
dprintf('############################# c_i classes ###############################');
c_class_info = suff_classify(max(fseq,[],1), fseq, HEL27_ccpt(:,1), options);
c_class_info.featfun = featfun;
c_class_info.time = toc;

if (clustercs)    
    c_class_info.cluster_info = class_info_cluster(c_class_info);
    
    ncl = max(c_class_info.cl(:,end));
    mincl = ceil(0.2*ncl);
    maxcl = ceil(0.8*ncl);
    c_class_info = class_info_applyclustering( c_class_info, [0.10,0.40], [mincl,maxcl], 50, 1e-3, 0.02 );
    cluster_info_plot(c_class_info.cluster_info, c_class_info.cluster_result);
    clc = c_class_info.cluster_result.cl;
else
    clc = c_class_info.cl(:,end);
end

%% cp classification
tic;
% find prediction for c_i+1                
dprintf('############################# c_i+1 classes ###############################');
cp_class_info = suff_classify(max(fseq,[],1), fseq, HEL27_ccpt(:,2), options);
cp_class_info.featfun = featfun;
cp_class_info.time = toc;

if (clustercs)
    cp_class_info.cluster_info = class_info_cluster(cp_class_info); 
    
    ncl = max(cp_class_info.cl(:,end));
    mincl = ceil(0.2*ncl);
    maxcl = ceil(0.8*ncl);
    cp_class_info = class_info_applyclustering( cp_class_info, [0.10,0.40], [mincl,maxcl], 50, 1e-3, 0.02 );
    cluster_info_plot(cp_class_info.cluster_info, cp_class_info.cluster_result);
    clcp = cp_class_info.cluster_result.cl;
else
    clcp = cp_class_info.cl(:,end);
end

%% final classification
tic;
dprintf('############################# final classes ###############################');
% generate pair featues from (seq(1:4), c_i, c_i+1)
triplepairfeats = feat_pairs([fragaaseq, clc, clcp]);

triples = getfreq([clc, clcp, hel27]);
figure;
hist(triples(:),50);
title 'triples histogram';

[trsort, perm] = sort([triples(:)]);
if (numtriples>0),
    ind = perm(end:-1:end-numtriples+1);
else
    ind = find(triples>mintriples);
end

s = sys10toX(ind,size(triples));
final_target_cl = zeros(size(hel27));
for c=1:size(s,1),
    ind = c_class_info.cl(:,end) == s(c,1) & cp_class_info.cl(:,end) == s(c,2) & hel27 == s(c,3);
    final_target_cl(ind) = c;
end

% classify 
options.new_classes = 0.015;
options.cat_ent_weight = 0.1;
if (max(final_target_cl)>30)
    options.conf_fn = @(x) x>0;
else
    if isfield(options,'conf_fn'),
        options = rmfield(options,'conf_fn');
    end
end
options.term_percentage = 0.025;
%options.savedata = 1; % needed for clustering
final_class_info = suff_classify(max(triplepairfeats,[],1), triplepairfeats, final_target_cl, options);
final_class_info.time = toc;

% hierarchical clustering of final classes
if clustert,    
    final_class_info.cluster_info = class_info_cluster(final_class_info);
    
    ncl = max(final_class_info.cl(:,end));
    mincl = ceil(0.15*ncl);
    maxcl = ceil(0.5*ncl);
    final_class_info = class_info_applyclustering( final_class_info, [0.10,0.50], [mincl,maxcl], 50, 1e-3, 0.015 );
    cluster_info_plot(final_class_info.cluster_info, final_class_info.cluster_result, hel27, final_class_info.cl(:,end));
    clf = final_class_info.cluster_result.cl;
else
    clf = final_class_info.cl(:,end);
end

%% cleanup
pr_comments = sprintf(input('post-run comments: ','s'));

new_deal_data.comments = comments;
new_deal_data.pr_comments = pr_comments;
new_deal_data.c = c_class_info;
new_deal_data.cp = cp_class_info;
new_deal_data.final = final_class_info;

hel27_quality(clf, hel27, 1);

% save diary
if (usediary)
    diary off;    
    fid = fopen('HEL27.dry','r');
    r = [];
    while 1,
        l = fgetl(fid);
        if ~ischar(l), 
            break;
        end
        s = sprintf('%s\n',l);
        r = [r s];
    end
    fclose(fid);
    rundiary = r;
    delete('HEL27.dry');
end
new_deal_data.diary = rundiary;

% save new_deal_data
global classdatadir;
save([classdatadir,filesep,'HEL27_gamma'],'new_deal_data');

% display results
HEL27_info('HEL27_gamma')
