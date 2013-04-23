% TESTSET Test classify on various data sets.
function [correct] = testset(setname, cvnfolds, doallfolds)


numdatasetdir = [fileparts(mfilename('fullpath')),filesep,'num_datasets'];
catdatasetdir = [fileparts(mfilename('fullpath')),filesep,'cat_datasets'];
mixeddatasetdir = [fileparts(mfilename('fullpath')),filesep,'mixed_datasets'];

if (nargin > 0)

    havenum = exist([numdatasetdir,filesep,setname,'.mat'],'file');
    havecat = exist([catdatasetdir,filesep,setname,'.mat'],'file');
    havemixed = exist([mixeddatasetdir,filesep,setname,'.mat'],'file');

    if (sum([havenum, havecat, havemixed]>0) > 1),
        dprintf('Ambiguous dataset name.');
        return;
    end
end

if nargin == 0 || ~any([havenum, havecat, havemixed]),
    dprintf('--------- Continuous data sets ---------');
    dir([numdatasetdir,filesep,'*.mat']);    
    dprintf('--------- Categorical data sets ---------');
    dir([catdatasetdir,filesep,'*.mat']);    
    dprintf('--------- Mixed data sets ---------');
    dir([mixeddatasetdir,filesep,'*.mat']);
    return;
end

if nargin < 2, cvnfolds = 10; end
if nargin < 3, doallfolds = 0; end

if (havenum), datasetdir =  numdatasetdir; end
if (havecat), datasetdir = catdatasetdir; end
if (havemixed), datasetdir = mixeddatasetdir; end

datafile = [datasetdir,filesep,setname,'.mat'];
load(datafile); % data, [datatype(in case of mixed)]; data(:,end) are classes
data = double(data);
[n, nattr] = size(data);

if (havenum), datatype = zeros(1,nattr-1); end
if (havecat), datatype = max(data(:,1:end-1),[],1); end

options = struct(   'max_iterations',1,...                    
                    'new_classes',0,...
                    'remove_small',0,...                    
                    'final_adjust',0,...  
                    'supervised',1,...
                    'autropy',3,...
                    'potential_options',suffclass.potential.defaultoptions,...
                    'potval_options',suffclass.potvals.defaultoptions);               
                
options.potval_options.use_entropy = 0;
options.potval_options.ent_weight = 1;
options.potval_options.handleequalpotentials = 0;
options.potval_options.use_potential_entropy = 1;
options.potential_options.debug_usematcov = 0;
options.potential_options.num_reg_weight = 0;

% check if we have premade test/train sets
trainname = regexprep(setname,'_test','_train');
if (~isempty(findstr('_test',setname)) && exist([datasetdir,filesep,trainname,'.mat'],'file'))  
    
    test_data = data(:,1:end-1);
    test_target = data(:,end);
    
    load([datasetdir,filesep,trainname]); % data, [datatype(in case of mixed)]
    train_data = data(:,1:end-1);
    train_target = data(:,end);
    testsize = length(test_target);
    
    % the first data loaded was the test data, but we want all categories,
    % which might not be present in the test set, so reget datatype
    if (havecat), datatype = max([train_data;test_data],[],1); end
    
	class_info = suffclass.classify(datatype,train_data,train_target,options);
    
    potvals = suffclass.potential_values(class_info.pot,class_info.options.potval_options);
    clpred = potvals.calcvals(test_data,class_info.potvals_priors);        
    correct = sum(test_target == clpred)/length(test_target);
    
    correct_train = sum(train_target == class_info.cl(:,end))/length(train_target);
    
    rndcl = randperm(testsize);
    rndcl = test_target(rndcl);
    confrnd = suffclass.utils.getfreq([test_target rndcl]);                        
    randoms = trace(confrnd)/testsize;
else
        
    cvobj = suffclass.test.crossval(n, cvnfolds);        
    
    for cvnum=1:cvnfolds,
        [testidx,trainidx] = cvobj.foldidxset(cvnum);

        train_data = data(trainidx,1:end-1);
        train_target = data(trainidx,end);

        test_data = data(testidx,1:end-1);
        test_target = data(testidx,end);

        testsize = length(test_target);
        
        class_info = suffclass.classify(datatype,train_data,train_target,options);

        % get percentage of correct train data + confus        
        correct_train(cvnum) = sum(train_target == class_info.cl(:,end))/length(train_target);
        
        % get percentage of correct test data + confus
        potvals = suffclass.potential_values(class_info.pot,class_info.options.potval_options);
        clpred = potvals.calcvals(test_data,class_info.potvals_priors);   
        correct(cvnum) = sum(test_target == clpred)/length(test_target); 
        
        % get percentage of random data
        rndcl = randperm(testsize);
        rndcl = test_target(rndcl);
        confrnd = suffclass.utils.getfreq([test_target rndcl]);                        
        randoms(cvnum) = trace(confrnd)/testsize;
        
        if (~doallfolds),     

            disp('training data confusion matrix:');
            confuspred = suffclass.utils.getfreq([train_target class_info.cl(:,end)]);
            suffclass.utils.tightmat(confuspred);
            q = sum(diag(confuspred))/sum(sum(confuspred));
            dprintf('diagonality: %.2f%%\n',q*100);

            disp('test data confusion matrix:');

            potvals = suffclass.potential_values(class_info.pot,class_info.options.potval_options);
            clpred = potvals.calcvals(test_data,class_info.potvals_priors);        
            confuspred = suffclass.utils.getfreq([test_target clpred]);
            suffclass.utils.tightmat(confuspred);
            q = sum(diag(confuspred))/sum(sum(confuspred));
            dprintf('diagonality: %.2f%%\n',q*100); 

            correct = trace(confuspred)/n;

            break;
        end
        
    end
end
    
disp('random class confusion matrix:');
%rndcl = fix(rand(testsize,1)*max(train_target)) + 1;
rndcl = randperm(testsize);
rndcl = test_target(rndcl);
confrnd = suffclass.utils.getfreq([test_target rndcl]);
disp(confrnd);
q = trace(confrnd)/sum(sum(confrnd));
dprintf('diagonality: %.2f%%',q*100);

cvnfolds = length(correct);
dprintf('CV with %i folds results:',cvnfolds);
dprintf('Mean+-std percentage correctly classified test data: %.2f+-%g',mean(correct)*100,std(correct*100));
dprintf('Mean+-std percentage correctly classified train data: %.2f+-%g',mean(correct_train)*100,std(correct_train*100));
dprintf('Mean+-std percentage randomly classified data: %.2f+-%g',mean(randoms)*100,std(randoms*100));