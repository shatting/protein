%% options
bdiagbounds = 1;
usediary = 1;
docluster = 0;
featfun = @feat_pairs;
%featfun = @(x) x;
options = struct(   'plots',1,...
                    'max_iterations',inf,...
                    'ask',0,...
                    'ask_timeout',5,...
                    'term_percentage',0.02,...
                    'new_classes',0.015,...
                    'remove_small',0.01,...                    
                    'stop_on_nonewclasses',0,...
                    'num_reg_weight',1,...
                    'num_ent_weight',1,...
                    'cat_prior_reg',0.01,...
                    'cat_cond_reg',0.05,...
                    'cat_ent_weight',1,...
                    'cat_no_lookup',1,...
                    'final_adjust',0,...
                    'savedata',0,...
                    'ent_weight',0.8,...
                    'handleequalpotentials',1,...
                    'supervised',1,...
                    'autropy',0);

%% preliminaries
clc;
close all;
diary off;
if (exist('HEFG.dry','file')), delete('HEFG.dry'); end;

seq = data.FeatureDB.db.getdata({'aa1','aa2','aa3','aa4'});
sclassifier = class.HEFGSClassifier(bdiagbounds);
hefg = sclassifier.classify(data.FeatureDB.db);

comments = sprintf(input('pre-run comments: ','s'));

if (usediary)                
    diary HEFG.dry
end
               
fseq = featfun(seq);
%fseq = fseq(1:1000,:);

tic;
% find prediction for hefg
dprintf('############################# hefg classes ###############################');

hefg_class_info = suffclass.classify(max(fseq,[],1), fseq, hefg, options);
hefg_class_info.featfun = featfun;
hefg_class_info.time = toc;

if (docluster)    
    options = hefg_class_info.options;
    options.clust_minclfact = 0.2;
    options.clust_maxclfact = 0.8;
    options.clust_relcostbds = [0.10, 0.40];
    options.clust_maxdiff = 50;
    options.clust_minprev = 1e-3;
    options.clust_minstep = 0.02;    
    hefg_class_info.options = options;
    hefg_class_info = class_info_cluster(hefg_class_info);   
end

new_deal_data = struct;
new_deal_data.comments = comments;
new_deal_data.hefg_info = hefg_class_info;
new_deal_data.diagbds = bdiagbounds;
new_deal_data.timedate = datestr(now,'dd.mm.yy, HH:MM:SS');

% assessment before comments
%HEFG_assessment

pr_comments = sprintf(input('post-run comments: ','s'));
new_deal_data.pr_comments = pr_comments;

% save diary
if (usediary)
    diary off;    
    fid = fopen('HEFG.dry','r');
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
    delete('HEFG.dry');
end
new_deal_data.diary = rundiary;

% save new_deal_data
k = 1;
global datadir;
fname = [datadir,filesep,'HEFG_gamma_',datestr(now,'yymmdd'),'.mat'];    
while exist(fname,'file')
    fname = [datadir,filesep,'HEFG_gamma_',datestr(now,'yymmdd'),'_',num2str(k),'.mat'];    
    k = k+1;    
end
save(fname,'new_deal_data');

% display results
%new_deal_info('new_deal_data')
new_deal_data
dprintf('saved to %s',fname);