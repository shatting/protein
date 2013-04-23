% opti.predictscript
%   frontend to opti.predict()
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
function [allpres,rrmsds,len] = predict_chain(optidata, hefg_tr, predict_options)

p_chain = data.GeomDB.ch(predict_options.chain_idx); % chain to predict
len = length(p_chain.seq);

dprintf('Running prediction on chain %s (%i in geomdb, %i aas)',...
    p_chain.name,predict_options.chain_idx, len);

% check if there are any aa > 20
if (any(p_chain.seq>20))
    error('chain has aas > 20');
else
    dprintf('ok, chain has no aas > 20.');
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

%opti.summarize_predict_result(allpres,rrmsds);