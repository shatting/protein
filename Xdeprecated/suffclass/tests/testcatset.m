% TESTCATSET Test suff_classify on various categorical data sets.
function testcatset(setname, nocv)

if nargin == 0 || ~exist(['cat_datasets',filesep,setname,'.mat'],'file'),
    dir(['cat_datasets',filesep,'*.mat']);
    return
end

if nargin == 1,
    nocv = 0;
end

load(['cat_datasets',filesep,setname]);

[n, nattr] = size(data);

if nocv,
    train_data = data(:,1:end-1);
    train_target = data(:,end);
    testsize = n;
    test_target = train_target;
else   
    cv = 10;
    testsize = fix(n/cv)
    trainsize = n - testsize
    perm = randperm(n);
    train_data = data(perm(1:trainsize),1:end-1);
    train_target = data(perm(1:trainsize),end);
    test_data = data(perm(trainsize+1:end),1:end-1);
    test_target = data(perm(trainsize+1:end),end);
end

options = struct(   'plots',0,...
                    'max_iterations',10,...
                    'ask',0,...
                    'ask_timeout',5,...
                    'term_percentage',0.02,...
                    'new_classes',0,...
                    'remove_small',0,...
                    'stop_on_nonewclasses',0,...
                    'num_reg_weight',1,...
                    'num_ent_weight',1,...
                    'cat_prior_reg',0.1,...
                    'cat_cond_reg',0.1,...
                    'cat_ent_weight',0,...
                    'cat_no_lookup',1,...
                    'final_adjust',0,...
                    'use_entropy',1,...
                    'handleequalpotentials',1,...
                    'savedata',0);
                
class_info = suff_classify(max(train_data,[],1),train_data,train_target, options);


if ~nocv,
    clpred = suff_data2v(test_data,class_info.pot,class_info.pot.catpot.ent,class_info.options);
    
    disp('test data confusion matrix:');
    confuspred = getfreq([test_target clpred]);
    disp(confuspred);
    q = sum(diag(confuspred))/sum(sum(confuspred));
    dprintf('diagonality: %.2f%%\n',q*100);
else
    disp('final confusion matrix:');
    disp(confus(:,:,end));
    q = trace(confus(:,:,end))/sum(sum(confus(:,:,end)));
    dprintf('diagonality: %.2f%%\n',q*100);
end
  
disp('random class confusion matrix:');
%rndcl = fix(rand(testsize,1)*max(train_target)) + 1;
rndcl = randperm(testsize);
rndcl = test_target(rndcl);
confrnd = getfreq([test_target rndcl]);
disp(confrnd);
q = trace(confrnd)/sum(sum(confrnd));
dprintf('diagonality: %.2f%%',q*100);
