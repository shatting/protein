% opti.predictscript
%   frontend to opti.predict(). dev version of predict_chain
%
%   with parallel computing toolbox:
%       run "matlabpool start x" where x is number of cores/workers

% [chs,idx] = data.GeomDB.db.getsizencell(10,15)

%  .options
%       .chain_idx      index of chain to predict in data.GeomDB
%       .HMM_numseqs    number of s-class sequences to initially try out
%       .useX0          whether to start optimization at actual fold
%       .usehydroconstr use hydro constraints
%       .usevmaxconstr  use V_max constraints
%       .useccconstr    use close contact constraints
%       .hefg_tr_id     id of HEFGTrainingRun to use

redo_optidata = 0; % always recalculate optidata?

predict_options = struct;

% ---  options concerning prediction
predict_options.chain_idx = 1139;
predict_options.HMM_numseqs = 100;
predict_options.useX0 = 0;
predict_options.usehydroconstr = 0;
predict_options.usevmaxconstr = 1;
predict_options.Vmaxpercentile = 0.9; % 0 = off
predict_options.useccconstr = 1;
% switch and factor in one
predict_options.useccconstr_useactualccmax = 1.3;
% use matlab covariance estimator for s-g combination covmodel
predict_options.use_sgmatcov = 1;

%'none' = knitro findiff, 'hand' = hand AD, 'intlab' = intlab overloaded AD
predict_options.gradient = 'hand'; 

% --- which potential to use
predict_options.hefg_tr_id = '130207_4';

% ------------- output options
predict_options.plotall = 0;
predict_options.savexhistory = 0;
predict_options.showsactualrmsdprogression = 0;
predict_options.getknitrotime = 0;

% ------------- performance options
predict_options.sactual = 1;     
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

p_chain = data.GeomDB.ch(predict_options.chain_idx); % chain to predict

dprintf('Testing prediction on protein %s (%i in geomdb, %i aas)',...
    p_chain.name,predict_options.chain_idx, length(p_chain.seq));

% check if there are any aa > 20
if (any(p_chain.seq>20))
    error('chain has aas > 20');
else
    dprintf('ok, chain has no aas > 20.');
end

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

% get correct geometry in terms of c and z
predict_options.X0 = [p_chain.geometry.c; p_chain.ztransform(hefg_tr.get_s_classifier())];

% CC constraint data
if (predict_options.useccconstr)
    if (predict_options.useccconstr_useactualccmax)    
        predict_options.ccmax = opti.constr_cc_findmax(p_chain)*...
            predict_options.useccconstr_useactualccmax;
    else        
        predict_options.ccmax = polyval(polycoeff,length(p_chain.seq));
    end
else
    predict_options.ccmax = 0;
end

% HYDRO constraint data
if (predict_options.usehydroconstr)    
    [ pot, chi ] = opti.constr_hydro_getpot(0);
    predict_options.hydro_chin = chi(length(p_chain.seq));
    predict_options.hydro_min = 0.7*opti.constr_hydro(p_chain.coords,...
        pot(p_chain.seq),chi(length(p_chain.seq)));
    predict_options.hydro_max = 1.3*opti.constr_hydro(p_chain.coords,...
        pot(p_chain.seq),chi(length(p_chain.seq)));
    predict_options.hydro_pot = pot;
end

% run prediction
[allpres, rrmsds, actualsprob, HMMseqs] =...
    opti.predict( p_chain, optidata, 'unused', predict_options );

opti.summarize_predict_result(allpres,rrmsds);