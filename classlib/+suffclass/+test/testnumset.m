% TESTCATSET Test suff_classify on various categorical data sets.
function [] = testnumset(setname, nocv)

datasetdir = [fileparts(mfilename('fullpath')),filesep,'num_datasets'];
if nargin == 0 || ~exist([datasetdir,filesep,setname,'.mat'],'file'),    
    dir([datasetdir,filesep,'*.mat']);
    return
end

if nargin == 1,
    nocv = 0;
end

load([datasetdir,filesep,setname]); % datatype, data
[n, nattr] = size(data);
trainname = regexprep(setname,'test','train');
if (~isempty(findstr('_test',setname)) && exist([datasetdir,filesep,trainname,'.mat'],'file'))  
    test_data = data(:,1:end-1);
    test_target = data(:,end);
    load([datasetdir,filesep,trainname]); % datatype, data
    train_data = data(:,1:end-1);
    train_target = data(:,end);
    testsize = length(test_target);
else
  
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
end

options = struct(   'max_iterations',1,...                    
                    'new_classes',0,...
                    'remove_small',0,...                    
                    'final_adjust',0,...  
                    'supervised',1,...
                    'autropy',0,...
                    'potential_options',suffclass.potential.defaultoptions,...
                    'potval_options',suffclass.potvals.defaultoptions);               
                
options.potval_options.use_entropy = 0;
options.potval_options.ent_weight = 1;
options.potval_options.handleequalpotentials = 1;
                
class_info = suffclass.classify(zeros(1,nattr-1),train_data,train_target,options);

if ~nocv,
    potvals = suffclass.potential_values(class_info.pot,class_info.options.potval_options);
    clpred = potvals.calcvals(test_data,class_info.potvals_priors);
    
    disp('training data confusion matrix:');
    confuspred = suffclass.utils.getfreq([train_target class_info.cl(:,end)]);
    suffclass.utils.tightmat(confuspred);
    q = sum(diag(confuspred))/sum(sum(confuspred));
    dprintf('diagonality: %.2f%%\n',q*100);
    
    disp('test data confusion matrix:');
    confuspred = suffclass.utils.getfreq([test_target clpred]);
    suffclass.utils.tightmat(confuspred);
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
confrnd = suffclass.utils.getfreq([test_target rndcl]);
disp(confrnd);
q = trace(confrnd)/sum(sum(confrnd));
dprintf('diagonality: %.2f%%',q*100);
