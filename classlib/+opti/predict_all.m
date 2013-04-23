redo_optidata = 0; % always recalculate optidata?

predict_options = struct;

% ---  options concerning prediction
predict_options.HMM_numseqs = 1;
predict_options.useX0 = 0;
predict_options.usehydroconstr = 0;
predict_options.usevmaxconstr = 0;
predict_options.Vmaxpercentile = 0.9; % 0 = off
predict_options.useccconstr = 0;
% switch and factor in one
predict_options.useccconstr_useactualccmax = 1.3;
% use matlab covariance estimator for s-g combination covmodel
predict_options.use_sgmatcov = 1;

%'none' = knitro findiff, 'hand' = hand AD, 'intlab' = intlab overloaded AD
predict_options.gradient = 'hand'; 

% --- which potential to use
%predict_options.hefg_tr_id = '130420';
predict_options.hefg_tr_id = '130207_4';

% ------------- output options
predict_options.plotall = 1;
predict_options.savexhistory = 0;
predict_options.showsactualrmsdprogression = 0;
predict_options.getknitrotime = 0;

% ------------- performance options
predict_options.sactual = 0;     
predict_options.fastobjective = 1; % not compatible with savexhistory
predict_options.priority = 'sn'; % process priority. see /library/priority

% --- options concerning optidata

% --- knitro options
predict_options.ktroptions=optimset('Display','off'); 
% final, off, iter-detailed, notify, notify-detailed,final-detailed

%predict_options.ktroptions=optimset(predict_options.ktroptions,'Algorithm','interior-point');
predict_options.ktroptions=optimset(predict_options.ktroptions,'Algorithm','active-set');
predict_options.ktroptions=optimset(predict_options.ktroptions,'MaxIter',200);

if (~strcmp(predict_options.gradient,'none'))
   % ktrlink wants both or none
   predict_options.ktroptions=optimset(predict_options.ktroptions,'GradObj','on','GradConstr','on');
else
   % finite differencing: central, forward
   predict_options.ktroptions=optimset(predict_options.ktroptions,'FinDiffType','forward'); 
end
%options.ktroptions=optimset(options.ktroptions,'PlotFcns',{@optimplotfirstorderopt}); 
% @optimplotx, @optimplotfunccount, @optimplotfval

% never, gradients in parallel
predict_options.ktroptions=optimset(predict_options.ktroptions,'UseParallel','always'); 

% set these to override default intlab behaviour
sparsegradient(1000);
feature jit on; % intlab turns these off
feature accel on;

% ------------- end of options

% check if optidata need to be rebuilt
if (~exist('current_hefg_tr_id','var') || ~exist('optidata','var') ||...
        (optidata.Vmaxpercentile ~= predict_options.Vmaxpercentile) || ...
         ~strcmp(predict_options.hefg_tr_id,current_hefg_tr_id) || redo_optidata)
    
    dprintf('Loading gamma training run "%s"',predict_options.hefg_tr_id);
    
    current_hefg_tr_id = predict_options.hefg_tr_id;
    
    hefg_tr = protein.HEFG.HEFGTrainingRun.load(current_hefg_tr_id);
    % TODO: remove the chain to predict from featuredb?
    optidata = opti.Optidata( data.FeatureDB.db, ...
        hefg_tr.get_s_classifier(), ...
        hefg_tr.get_g_classifier(), ...
        predict_options.hefg_tr_id, predict_options );
    
else
     dprintf('Using current gamma training run "%s"', predict_options.hefg_tr_id);
end

allpres = cell(1334,1);
rrmsds = zeros(1334,predict_options.HMM_numseqs);
len = zeros(1334,1);
for chnum = 1:1334
    try
        predict_options.chain_idx = chnum;
        [~,rrmsds(chnum,:), len(chnum)] = opti.predict_chain(optidata, hefg_tr, predict_options);
    end
end
