% SUFF_CLASSIFY Partially or fully supervised iterative classification of mixed data.
%
% cl_info = suff_classify(type, data, target, options)
%
% target labels are known classes if positive, unknown if 0
%
% The predicted classes may have different labels than initially.
% Predict original target labels by label=argmin(cost*confus(:,cl))
% where cost(s,t) is the cost of deciding for label s given target t
%
% INPUT:
% type(i)               type of feature i (see suffstat.m)
% data(l,i)             i-th feature of datapoint l
% cl(l,1)               initial class assignment
%                       unknown if 0
% options               .plot                   should plots be done (0)
%                       .max_iterations         maximum iterations (inf)
%                       .ask                    ask after each iteration(0)
%                       .ask_timeout            timeout for ask(5)
%                       .term_percentage        terminate if
%                                               relabeled<this (0.05)
%                       .new_classes            what percentage of
%                                               total data points should be
%                                               exceeded in iteration
%                                               confusion matrix to form
%                                               new class (0.025)
%                       .remove_small           
%                       .stop_on_nonewclasses   stop if no new classes are
%                                               formed (0)
%                       .num_reg_weight         weight of total potential
%                                               in regularization of
%                                               numerical class potentials
%                                               (1)
%                       .num_ent_weight         weight of entropy term in 
%                                               numerical part (1)
%                       .cat_prior_reg          number of dummy categorical             
%                                               class points (0.01)
%                       .cat_cond_reg           number of dummy categorical
%                                               conditional points
%                                               (realmin)
%                       .cat_ent_weight         weight of categorical
%                                               priors (0.5). 0 for uniform
%                       .cat_no_lookup          use cat lookup? (0)
%                       .final_adjust           adjust potential after
%                                               stopping to reflect final
%                                               class changes? (0)
%                       .savedata               data saved to output? (0)
%                       .supervised             compare to initial classes
%                                               in each iteration (0)
%                                               .

% 
% OUTPUT:
% cl_info               classification information struct with fields
%   .options            options that were used
%   .pot:               see classpot.m
%   .freqreg:           see classpot.m
%   .potmin:            see classpot.m
%   .prob:              see classpot.m    (savedata=1)
%   .cl(l,it)           class assignment of l-th data point in iteration it
%   .confus             start-end confusion matrix
%   .freq               non-regularized frequencies (in contrast to
%                       .freqreg)
%   .data:              n x dim data used (savedata=1)
%   .suff:              suff that was used to generate .pot, (savedata=1)
%   .type:              1 x d data type (<0,==0,>0), see suffstat.m

function [cl_info]=suff_classify(type, data, cl, options)

dprintf('\n\n-------------- start of suff_classify ------------');
% options
if nargin<4, options = struct; end
if ~isfield(options,'plots'),           options.plots=0; end;
if ~isfield(options,'max_iterations'),  options.max_iterations=inf; end;
if ~isfield(options,'ask'),             options.ask=0; end;
if ~isfield(options,'ask_timeout'),     options.ask_timeout=5; end; % use 0 for no timeout
if ~isfield(options,'term_percentage'), options.term_percentage = 0.05; end;

if ~isfield(options,'new_classes'),     options.new_classes = 0.025; end;
if ~isfield(options,'remove_small'),    options.remove_small = 0.025; end;
if ~isfield(options,'stop_on_nonewclasses'),     options.stop_on_nonewclasses = 0; end;

if ~isfield(options,'num_reg_weight'),  options.num_reg_weight = 1; end;
if ~isfield(options,'num_ent_weight'),  options.num_ent_weight = 1; end;

if ~isfield(options,'cat_prior_reg'),   options.cat_prior_reg = 0.01; end;
if ~isfield(options,'cat_cond_reg'),    options.cat_cond_reg = realmin; end;
if ~isfield(options,'cat_ent_weight'),  options.cat_ent_weight = 0.5; end;

if ~isfield(options,'cat_no_lookup'),   options.cat_no_lookup = 0; end;
if ~isfield(options,'final_adjust'),    options.final_adjust= 0; end;
if ~isfield(options,'savedata'),        options.savedata= 0; end;
if ~isfield(options,'autropy'),         options.autropy= 0; end;
if ~isfield(options,'supervised'),      options.supervised= 0; end;

if (options.autropy)
    if (exist('autropy','dir') == 7)
        addpath('autropy');
    else
        warning('autropy chosen by option, but package directory not found. not using autropy.');
        options.autropy = 0;
    end
end

dprintf('options:');
disp(options);

[nd,dim] = size(data);    % nd items of length dim
dprintf('%i data points of dimension %i.',nd, dim);
if nd==0, error('empty data set'); end;

% initial class assignment
ncl = max(cl);
nnocl = sum(cl==0);
dprintf('%i target classes.',ncl);
if nnocl>0,
    dprintf('%i unlabelled items.',nnocl);
else
    dprintf('no unlabelled items.');
end

freq=getfreq(cl);

dprintf('feature types: [%s]',sprintf('%i, ',type)); 

suff = suffstat(ncl,type);

%check for 0 cat
m = min(data(:,suff.cat),[],1);
if (any(m<=0)),
    error('categorical features must be > 0.');
end

nreassignments(1:10)=inf;
iteration=0;

cl_idx_lab_1 = cl(:,1) > 0; 
cl_idx_unlab_1 = cl(:,1) == 0;

while 1,
    tic;
    iteration=iteration+1;  
    ifreq = getfreq(cl(:,iteration));    
    incl = max(cl(:,iteration));
    cl_idx_unlab_i = cl(:,iteration)==0;
    cl_idx_lab_i = cl(:,iteration)>0;

    dprintf('____________________________________________________________');
    dprintf('ITERATION %i\n', iteration);
    dprintf('%i classes, frequencies:', incl);  
    tightmat(ifreq);
    if (any(cl_idx_unlab_i))
        dprintf('%i unlabelled items.', sum(cl_idx_unlab_i));
    else
        dprintf('no unlabelled items.');
    end

    % sort classes by frequency
    if 0,
      [nfsort,perm]=sort(ifreq);
      invp(perm)=[ncl:-1:1]';
      cl(:,iteration)=invp(cl(:,iteration));

      % merge tiny classes 
      if iteration>1,   
        large=(sum(ifreq.^2>=nd));
        if large<ncl, 
          havemerged = 1;
          ncl=1+sum(ifreq.^2>=nd);
          cl(cl(:,iteration)>ncl,iteration)=ncl;
        end;
      end
    end;

    % build sufficient statistic  
    %   if (iteration > 1 & 0),
    %       dprintf('adjusting sufficient statistics');
    %       if (oldncl < ncl)
    %           suff = suffstat(suff,ncl);
    %       end
    %       suff = suffstat(suff,data,nd,cl(:,iteration),cl(:,iteration-1));
    %   else      
    % removed adjustment because it seems slower and not sure if its
    % correct
	dprintf('building sufficient statistics and potential..');
    suff = suffstat(incl,type);
    suff = suffstat(suff,data,size(data,1),cl(:,iteration)); % ignores cl==0
    %  end

    pot = suff_pot( suff, options );

    ifreqrel = ifreq/sum(ifreq); % /sum() instead of /nd because of possible empty classes  
    if (options.autropy),
      if (~min(cl(:,end)))
          error('partially supervised autropy not yet tested');
      end
      Vgd = vgammadelta( pot, data, cl(:,iteration) );

      [ ent, delta ] = getentropies( Vgd )
      options.use_entropy = 1;

      data2vpriors = ent;
    else      
      options.use_entropy = 0;
      data2vpriors = ifreqrel;      
    end

    dprintf('finding preliminary classes under potential.\n');
    [ cl_potential, potmin, potvals] = suff_data2v( data, pot, data2vpriors, options );

    % show initial-preliminary confusion matrix 
    showconfus(cl(cl_idx_lab_1,1),cl_potential(cl_idx_lab_1),'initial-(current potential) assignment confusion matrix (excluding initially unlabelled), rows index=c(1), column index=classes under pot');

    % show initially unlabelled preliminary assignment
    if sum(cl_idx_unlab_1) > 0,
     dprintf('initially unlabelled potential assignment');
     tightmat(getfreq([cl(cl_idx_unlab_1,1)+1 cl_potential(cl_idx_unlab_1)]));
    end

    if (iteration>1)
        % show previous iteration to preliminary assignment confusion matrix                    
        str = sprintf('before to after potential classification confusion matrix (excluding previously unlabelled), row index=c(%i), column index=classes under pot',iteration);
        showconfus(cl(cl_idx_lab_i,iteration),cl_potential(cl_idx_lab_i),str);

        % previously unlabelled items    
        if (any(cl_idx_unlab_i)),        
            dprintf('previously unlabelled items potential assignment:');
            tightmat(getfreq([cl(cl_idx_unlab_i,iteration)+1 cl_potential(cl_idx_unlab_i)]));    
        end        
    end
    dprintf('before and after potential class distributions (1st row before, 2nd after, 3rd after only previously unlabelled)');
    % we might have lost our last class
    f2 = getfreq(cl_potential);
    if length(f2) < max(cl(:,iteration))
        f2(max(cl(:,iteration))) = 0;
    end
    % unlabelled items might not have been assigned to the last class
    f3 = getfreq(cl_potential(cl_idx_unlab_i));
    if length(f3) < max(cl(:,iteration))
        f3(max(cl(:,iteration))) = 0;
    end
    tightmat([getfreq(cl(:,iteration));f2;f3] );    
        
    % at this stage, we do not have unlabelled items, but we may have empty
    % classes -> remove them.

    ifreq_pot = getfreq(cl_potential);
    ncl_pot = length(ifreq_pot); % number of classes hasnt changed yet, though..
    ifreq_pot_rel = ifreq_pot/nd; % /nd because we want frequencies relative to whole data set, ie including empty classes
    idx_ifreq_classestokeep = ifreq_pot_rel >= options.remove_small; 
    nclassesremove = ncl_pot - length(idx_ifreq_classestokeep);
    if (nclassesremove)      
      % expand clremove
      clremove(iteration,ncl_pot) = 0;            
      % find assignment
      clremove(iteration,idx_ifreq_classestokeep) = 1:length(idx_ifreq_classestokeep);

      dprintf('removing %i class(es) from potential assignment (first row = label before, second row = label after):',nclassesremove);
      tightmat([1:ncl_pot; clremove(iteration,:)]);            

      % apply assignment
      cl_afterremoval = clremove(iteration,cl_potential)';
      dprintf('%i items are now unlabelled.',nd - sum(cl_afterremoval~=0));
    else
        cl_afterremoval = cl_potential;
    end

    % cl_afterremoval does not contain empty classes 
    if (any(freqs(cl_afterremoval)==0)),
      error('there still are empty classes after removal!');
    end  
    % ... but may have unlabelled items (if we removed a class)    
    cl_afternew = cl_afterremoval;

    % find new classes. do not touch unlabelled.
    if options.new_classes,      
      dprintf('finding new classes..'); %TODO: make new class only if target potential class from which items are take away would not collapse
      
      cl_idx_unlab_afterremoval = cl_afterremoval == 0;         

      % confusion matrix used for determining new classes
      if (options.supervised || iteration == 1)
        previtidx = 1;
        cl_idx_new_classes = ~cl_idx_unlab_afterremoval & cl_idx_lab_1; % items that qualify for formation of new classes
      else
        previtidx = iteration;
        cl_idx_new_classes = ~cl_idx_unlab_afterremoval & cl_idx_lab_i; % items that qualify for formation of new classes
      end

      newclassconf = getfreq([cl(cl_idx_new_classes,previtidx) cl_afterremoval(cl_idx_new_classes)]);

      newclassconf_rel = newclassconf/nd;
      idx_newclassconf = (newclassconf_rel > options.new_classes) & (1-eye(size(newclassconf_rel))); % we dont want correctly classifieds as new classes                  
      nnewcl = sum(idx_newclassconf(:));      

      dprintf('entries with more than %i (%.2f%%) items get new class label.',floor(nd*options.new_classes),options.new_classes*100);

      if (nnewcl)    
          % display new classes in confusion style
          dprintf('%i new classes found [row=c(%i), column=c_pot, value=c(%i)]',nnewcl,previtidx,iteration+1);          
          newclassesmat = zeros(size(newclassconf_rel));
          newclassesmat(idx_newclassconf) = [1:nnewcl] + incl;
          tightmat(newclassesmat);                    

          % get new assignment
          [clprev, clpot] = find(idx_newclassconf);   
          
          % display new classes in table style
          dprintf('        c(%i)    c(%i)     c_pot     n', iteration+1,previtidx);            
          tightmat([[1:nnewcl]'+incl clprev clpot newclassconf(idx_newclassconf(:))],'%10.0f','%2i');

          % assign new labels                
          for i=1:nnewcl,
              ind = cl(:,previtidx) == clprev(i) & cl_afterremoval == clpot(i);
              cl_afternew(ind) = incl + i;
          end
      else
          dprintf('no new classes found.');
      end

      str = sprintf('initial-current iteration confusion matrix (excluding unlabelled), rows index=c(1), column index=c(%i)',iteration+1);
      showconfus(cl(cl_idx_new_classes,1),cl_afternew(cl_idx_new_classes),str);

      if (iteration>1)
         str = sprintf('previous-current iteration confusion matrix (excluding unlabelled), rows index=c(%i), column index=c(%i)',iteration,iteration+1);
         showconfus(cl(cl_idx_new_classes,iteration),cl_afternew(cl_idx_new_classes),str);
      end

    else
      nnewcl = 0;      
    end

    % done changing classes -> assign  
    cl(:,iteration+1) = cl_afternew;  
    ifreq = freqs(cl_afternew)'; % does not include unlabelled
    freq(iteration+1,1:length(ifreq)) = ifreq;  

    % stopping test
    itdur = toc;
    ndif=sum(cl(:,iteration)~=cl(:,iteration+1));    
    nunlab = sum(cl(:,iteration+1) == 0);

    dprintf('%i of %i (%.2f%%) items have new class after this iteration.',ndif,nd,ndif/nd*100);
    dprintf('%i of %i (%.2f%%) unlabelled items after this iteration.',nunlab,nd,nunlab/nd*100);    
    dprintf('%.2f s spent in iteration, an average of %.1f datapoints/s',itdur,nd/itdur);   
    dprintf('number of reassignments in last %i iterations:',length(nreassignments));
    tightmat([nreassignments(2:end) ndif]);

    % stopping tests
    bstop = 0;    
    if (ndif == 0), dprintf('no change in iteration.'); bstop = 1; end    
    if (ndif>=max(nreassignments)), disp('ndif>=max(nreassignments)'); bstop = 1; end    
    if (ndif < options.term_percentage*nd), dprintf('number of relabeled items < %f%%',options.term_percentage*100); bstop = 1; end
    if (iteration >= options.max_iterations), disp('iteration limit reached'); bstop = 1; end
    if (options.stop_on_nonewclasses && nnewcl > 0), disp('no new classes'); bstop = 1; end      
    if (options.ask && ~ask('continue?',options.ask_timeout)); bstop = 1; end      
    
    if bstop, 
        disp('stopping'); 
        break; 
    end    
    nreassignments=[nreassignments(2:length(nreassignments)),ndif];
end;

% adjust pot
%if ((ndif>0 && options.final_adjust) || options.remove_small && ~min(cl(:,end))) %% if done, the classes dont conform with pot!
    dprintf('calculating final potential and classes.');
    nclend = max(cl(:,end));
    suff = suffstat(nclend,type);
    suff = suffstat(suff,data,nd,cl(:,end));
    pot = suff_pot( suff, options );
    pr = freqs(cl(:,end))/sum(freqs(cl(:,end))); %TODO: or autropy
    [ cl(:,end+1), potmin, potvals] = suff_data2v( data, pot, pr, options );
    
    ifreq = freqs(cl(:,end))';
    freq(end+1,1:length(ifreq)) = ifreq;      
%end

dprintf('\ninitial & final class distributions:');
disp([freq(1,:);freq(end,:)]);

% assemble cl_info
cl_info.pot = pot;
if (max(cl) <= 255)
    cl_info.cl = uint8(cl);
else if max(cl) <= 255*255
    cl_info.cl = uint16(cl);        
    end
end
if (options.autropy)
    cl_info.ent = ent;
else
    cl_info.pr =pr;
end

cl_info.freq = freq;
if (options.savedata)    
    cl_info.data = data;
end
cl_info.suff = suff;

cl_info.options = options;

cl_info.readme =  [
'        pot: see suff_pot.m                '
'         cl: [ndata X niterations]         '
'       freq: [niterations X nclasses]      '
'       data: n x dim data used (savedata=1)'
'       suff: suff used for .pot            '
'     pr/ent: pr/ent used for suff_data2v() '
'       type: 1 x d data type (<0,==0,>0)   '
'    options: options used                  '
'     readme: this readme                   '];

dprintf('--------------end of suff_classify ------------\n\n');


    function [confus, diagp] = showconfus(cl1, cl2,dprintstr)
      confus = getfreq([cl1 cl2]);  
      diagp = sum(diag(confus))/sum(confus(:));

      dprintf(dprintstr);
      dprintf('diag=%.2f%%',diagp*100);
      
      if size(confus,2) < incl,
            confus(1,incl) = 0;
      end      
      
      tightmat([[sum(confus(:)) , -1 ,sum(confus,1)] ; zeros(1,size(confus,2)+2)-1 ; [sum(confus,2), zeros(size(confus,1),1)-1,confus]]);

    end

end