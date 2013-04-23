% suffclass.classify
%   Partially or fully supervised iterative classification of mixed data.
%
% cl_info = classify(type, data, target, options)
%
% target labels are known classes if positive, unknown if 0
%
% The predicted classes may have different labels than initially.
% Predict original target labels by label=argmin(cost*confus(:,cl))
% where cost(s,t) is the cost of deciding for label s given target t
%
% INPUT:
%   type(i)               type of feature i (see suffclass.suffstat)
%   data(l,i)             i-th feature of datapoint l
%   cl(l,1)               initial class assignment, unknown if 0
%
%   options
%       .plot               should plots be done (0)
%       .max_iterations     maximum iterations (inf)
%       .ask_timeout        timeout for ask(5). 0 to disable
%       .term_percentage    terminate if relabeled<this (0.05)
%       .new_classes        what percentage of total data points should be
%                           exceeded in iteration confusion matrix to form
%                           new class (0.025)
%       .remove_small       what percentage of total data points should be
%                           exceeded in iteration confusion matrix for the
%                           corresponding class to be removed (0.025)
%       .stop_on_nonewclasses   stop if no new classes are formed (0)
%       .final_adjust       adjust potential after stopping to reflect final
%                           class changes? (0)
%       .savedata           data(:,:) saved to output? (0)
%       .supervised         compare to initial classes in each iteration to
%                           determine new classes (0).
%                           leaving this at 0 will use last iteration's 
%                           classes
%       .autropy            use entropy optimization?
%       .potential_options  options to use for suffclass.potential
%       .potval_options     options to use for suffclass.potvals    
% 
% OUTPUT:
% cl_info               classification information struct with fields
%   .options            options that were used
%   .pot:               see suffclass.potential
%   .freqreg:           see suffclass.potential
%   .potmin:            see suffclass.potential
%   .prob:              see suffclass.potential    (savedata=1)
%   .cl(l,it)           class assignment of l-th data point in iteration it
%   .confus             start-end confusion matrix
%   .freq               non-regularized frequencies (in contrast to
%                       .freqreg)
%   .data:              n x dim data used (savedata=1)
%   .suff:              suff that was used to generate .pot, (savedata=1)
%   .type:              1 x d data type (<0,==0,>0), see suffstat.m

function [cl_info]=classify(type, data, cl, options)

% import suffclass.*;
% import suffclass.utils.*;
tic;

dprintf('\n\n-------------- start of suff_classify ------------');
if nargin<4, options = struct; end
if ~isfield(options,'plots'),           options.plots=0; end;
if ~isfield(options,'max_iterations'),  options.max_iterations=inf; end;
if ~isfield(options,'ask_timeout'),     options.ask_timeout=0; end;
if ~isfield(options,'term_percentage'), options.term_percentage = 0.05; end;
if ~isfield(options,'new_classes'),     options.new_classes = 0.025; end;
if ~isfield(options,'remove_small'),    options.remove_small = 0.025; end;
if ~isfield(options,'stop_on_nonewclasses'),  
                                    options.stop_on_nonewclasses = 0; end;
if ~isfield(options,'final_adjust'),    options.final_adjust= 0; end;
if ~isfield(options,'savedata'),        options.savedata= 0; end;
if ~isfield(options,'autropy'),         options.autropy= 0; end;
if ~isfield(options,'supervised'),      options.supervised= 0; end;

if ~isfield(options,'potential_options'),  options.potential_options=struct; end;
options.potential_options = ... 
    suffclass.potential.defaultoptions(options.potential_options);

if ~isfield(options,'potval_options'),  options.potval_options=struct; end;
options.potval_options = ...
    suffclass.potential_values.defaultoptions(options.potval_options);

% TODO: the following should go somewhere else
% if (options.autropy > 0 && options.potval_options.use_potential_entropy), 
%     disp('autropy selected, overriding potval_options.use_potential_entropy to false');
%     options.potval_options.use_potential_entropy = 0; 
% end

% check data
[nd,dim] = size(data);
if nd==0, error('empty data set'); end;

% display options
dprintf('--- options:');
disp(options);

dprintf('\n\n--- %i data points, %i features', nd, dim);
dprintf('--- feature types: [%s]',sprintf('%i, ',type));

% initial class assignment
% TODO: unique classes, remove empty ones
ncl = max(cl); % number of classe
nnocl = sum(cl==0); % number of unassigned
dprintf('--- %i initial classes, %i unlabelled items.',ncl,nnocl);

freq = suffclass.utils.getfreq(cl); % class frequencies
suffclass.display.tightmat(freq);
suffclass.display.tightmat(round(freq/sum(freq)*100));

nreassignments(1:10)=inf;
it=0; % iteration counter

cl_idx_lab_1 = cl(:,1) > 0; % index of initially labelled items
cl_idx_unlab_1 = cl(:,1) == 0; % index of initially unlabelled items

% main loop
while 1,
    tic;
    it=it+1;
    freq_i = suffclass.utils.getfreq(cl(:,it));  % class frequencies in current iteration
    ncl_i = max(cl(:,it));       % number of classes in current iteration
    cl_idx_unlab_i = cl(:,it)==0;    
    cl_idx_lab_i = cl(:,it)>0;
    
    dprintf('_____________ ITERATION %i ___________________________\n', it);
    dprintf('%i classes, frequencies:', ncl_i);  
    suffclass.display.tightmat(freq_i);
    if (sum(cl_idx_unlab_i)>0), dprintf('%i unlabelled items', sum(cl_idx_unlab_i)); end;

    % sort classes by frequency
    if 0,
      [nfsort,perm]=sort(freq_i);
      invp(perm)=[ncl:-1:1]';
      cl(:,it)=invp(cl(:,it));

      % merge tiny classes 
      if it>1,   
        large=(sum(freq_i.^2>=nd));
        if large<ncl, 
          havemerged = 1;
          ncl=1+sum(freq_i.^2>=nd);
          cl(cl(:,it)>ncl,it)=ncl;
        end;
      end
    end;
    
    % TODO: modify statistic instead of building it new. removed that part
    % because it was slower and not entirely correct
    dprintf('calculating suff statistics..');
    suff = suffclass.suffstat(ncl_i,type);    % create suffstat object
    suff.adddata(data,cl(:,it));     % add data with current class assigment
    
    dprintf('calculating potential..'); % create potential object    
    if (options.potential_options.debug_usematcov)
        pot = suffclass.potential( suff, options.potential_options, data, cl(:,it) );
    else
        pot = suffclass.potential( suff, options.potential_options);
    end
    % relative requencies. /sum() instead of /nd because of possible empty 
    % classes
    freq_i_rel = freq_i/sum(freq_i);
    
    % TODO: following should go into @potential_values
    if (options.autropy),
      if (any(cl_idx_unlab_i))
          error('partially supervised autropy not yet tested');
      end
      
      entropyobj = suffclass.entropyopt(pot, data, cl(:,it), options.potval_options);      
      [ data2vpriors, delta ] = entropyobj.getentropies(options.autropy)
           
      if (all(delta) > 0)                    
        dprintf('Entropy optimization found optimal separation.');        
      else
        dprintf('No optimal separation achieveable with entropy optimization.');                
      end
    else      
      data2vpriors = freq_i_rel;      
    end

    %dprintf('finding preliminary classes under potential\n');    
    potvals = suffclass.potential_values( pot, options.potval_options );
    cl_potential = potvals.calcvals( data, data2vpriors, ~options.autropy );
    
    % show initial-current confusion matrix 
    suffclass.display.showconfus(cl(cl_idx_lab_1,1), ...
        cl_potential(cl_idx_lab_1),ncl_i,...
        'initial-(current potential) assignment confusion matrix (excluding initially unlabelled), rows index=c(1), column index=classes under pot');

    % show current assignment of initially unlabelled items
    if sum(cl_idx_unlab_1) > 0,
     dprintf('assignment of initially unlabelled items');
     suffclass.display.tightmat(suffclass.utils.getfreq([cl(cl_idx_unlab_1,1)+1 ...
         cl_potential(cl_idx_unlab_1)]));
    end

    if (it>1)
        % show previous it to preliminary assignment confusion matrix                            
        suffclass.display.showconfus(cl(cl_idx_lab_i,it),...
            cl_potential(cl_idx_lab_i),ncl_i,...
            sprintf('before to after potential classification confusion matrix (excluding previously unlabelled), row index=c(%i), column index=classes under pot',it));

        % previously unlabelled items    
        if (any(cl_idx_unlab_i)),        
            dprintf('previously unlabelled items potential assignment:');
            suffclass.display.tightmat(suffclass.utils.getfreq([cl(cl_idx_unlab_i,it)+1 cl_potential(cl_idx_unlab_i)]));    
        end        
    end
    dprintf('before and after potential class distributions (1st row before, 2nd after, 3rd after only previously unlabelled)');
    % we might have lost our last class
    f2 = suffclass.utils.getfreq(cl_potential);
    if length(f2) < max(cl(:,it))
        f2(max(cl(:,it))) = 0;
    end
    % unlabelled items might not have been assigned to the last class
    f3 = suffclass.utils.getfreq(cl_potential(cl_idx_unlab_i));
    if length(f3) < max(cl(:,it))
        f3(max(cl(:,it))) = 0;
    end
    suffclass.display.tightmat([suffclass.utils.getfreq(cl(:,it));f2;f3] );    
        
    % at this stage, we do not have unlabelled items, but we may have empty
    % classes -> remove them.
    ifreq_pot = suffclass.utils.getfreq(cl_potential);
    ncl_pot = length(ifreq_pot); % number of classes hasnt changed yet
    
    % /nd because we want frequencies relative to whole data set, ie
    % including empty classes
    ifreq_pot_rel = ifreq_pot/nd; 
    idx_ifreq_classestokeep = ifreq_pot_rel >= options.remove_small; 
    nclassesremove = ncl_pot - length(idx_ifreq_classestokeep);
    if (nclassesremove)      
      % expand clremove
      clremove(it,ncl_pot) = 0;            
      % find assignment
      clremove(it,idx_ifreq_classestokeep) = ...
          1:length(idx_ifreq_classestokeep);

      dprintf(['removing %i class(es) from potential',...
        '(first row = label before, second row = label after):'],...
            nclassesremove);
      suffclass.display.tightmat([1:ncl_pot; clremove(it,:)]);            

      % apply assignment
      cl_afterremoval = clremove(it,cl_potential)';
      dprintf('%i items are now unlabelled.',nd - sum(cl_afterremoval~=0));
    else
        cl_afterremoval = cl_potential;
    end

    % cl_afterremoval should not contain empty classes 
    if (any(suffclass.utils.getfreq(cl_afterremoval)==0)),
      warning('there still are empty classes after removal!');
    end  
    % ... but may have unlabelled items (if we removed a class)    
    cl_afternew = cl_afterremoval;

    % find new classes. do not touch unlabelled.
    if options.new_classes,      
      %TODO: make new class only if target potential class from which items
      %are taken away would not collapse
      dprintf('finding new classes..'); 
      
      cl_idx_unlab_afterremoval = cl_afterremoval == 0;         

      % confusion matrix used for determining new classes
      if (options.supervised || it == 1)
        previtidx = 1;
        % items that qualify for formation of new classes
        cl_idx_new_classes = ~cl_idx_unlab_afterremoval & cl_idx_lab_1; 
      else
        previtidx = it;
        % items that qualify for formation of new classes
        cl_idx_new_classes = ~cl_idx_unlab_afterremoval & cl_idx_lab_i; 
      end

      newclassconf = suffclass.utils.getfreq([cl(cl_idx_new_classes,previtidx)...
          cl_afterremoval(cl_idx_new_classes)]);

      newclassconf_rel = newclassconf/nd;
      % we dont want correctly classifieds as new classes   
      idx_newclassconf = (newclassconf_rel > options.new_classes) & ...
          (1-eye(size(newclassconf_rel)));               
      nnewcl = sum(idx_newclassconf(:));      

      dprintf('entries with more than %i (%.2f%%) items get new class label.',...
          floor(nd*options.new_classes),options.new_classes*100);

      if (nnewcl)    
          % display new classes in confusion style
          dprintf('%i new classes found [row=c(%i), column=c_pot, value=c(%i)]',...
              nnewcl,previtidx,it+1);          
          newclassesmat = zeros(size(newclassconf_rel));
          newclassesmat(idx_newclassconf) = [1:nnewcl] + ncl_i;
          suffclass.display.tightmat(newclassesmat);                    

          % get new assignment
          [clprev, clpot] = find(idx_newclassconf);   
          
          % display new classes in table style
          dprintf('        c(%i)    c(%i)     c_pot     n', it+1,previtidx);            
          suffclass.display.tightmat([[1:nnewcl]'+ncl_i clprev clpot ...
              newclassconf(idx_newclassconf(:))],'%10.0f','%2i');

          % assign new labels                
          for j=1:nnewcl,
              ind = cl(:,previtidx) == clprev(j) & cl_afterremoval == clpot(j);
              cl_afternew(ind) = ncl_i + j;
          end
      else
          dprintf('no new classes found.');
      end

      str = sprintf('initial-current confusion matrix (excluding unlabelled), rows index=c(1), column index=c(%i)',it+1);
      suffclass.display.showconfus(cl(cl_idx_new_classes,1),cl_afternew(cl_idx_new_classes),ncl_i,str);

      if (it>1)
         str = sprintf('previous-current confusion matrix (excluding unlabelled), rows index=c(%i), column index=c(%i)',it,it+1);
         suffclass.display.showconfus(cl(cl_idx_new_classes,it),cl_afternew(cl_idx_new_classes),ncl_i,str);
      end

    else
      nnewcl = 0;      
    end

    % done changing classes -> assign  
    cl(:,it+1) = cl_afternew;  
    freq_i = suffclass.utils.getfreq(cl_afternew)'; % does not include unlabelled
    freq(it+1,1:length(freq_i)) = freq_i;  

    % stopping test
    itdur = toc;
    ndif=sum(cl(:,it)~=cl(:,it+1));    
    nunlab = sum(cl(:,it+1) == 0);

    nreass_new = [nreassignments(2:end) ndif];
    nreass_new(nreass_new==Inf)=[];
    dprintf('%i total / %i new classes',max(cl(:,it+1)),nnewcl)
    dprintf('%i of %i (%.2f%%) items have new class after this iteration (term: < %.2f%%)',...
        ndif,nd,ndif/nd*100,options.term_percentage*100);
    if (nunlab>0), dprintf('%i of %i (%.2f%%) unlabelled items after this iteration.',...
            nunlab,nd,nunlab/nd*100); end
    dprintf('number of reassignments in last %i iterations:',length(nreassignments));
    suffclass.display.tightmat(nreass_new);
    dprintf('percentage of reassignments in last %i iterations:',length(nreassignments));
    suffclass.display.tightmat(round(nreass_new/nd*100));
    dprintf('%.2f s spent in iteration, an average of %.1f datapoints/s',itdur,nd/itdur);   

    % stopping tests
    b_stop = 0;    
    if (ndif == 0), dprintf('no reassignments in iteration.'); b_stop = 1; end    
    if (ndif>=max(nreassignments)), disp('ndif>=max(nreassignments)'); b_stop = 1; end    
    if (ndif < options.term_percentage*nd), 
        dprintf('number of reassigned items < %f%%',options.term_percentage*100); b_stop = 1; 
    end
    if (it >= options.max_iterations), disp('iteration limit reached'); b_stop = 1; end
    if (options.stop_on_nonewclasses && nnewcl == 0), disp('no new classes'); b_stop = 1; end      
    if (~b_stop && options.ask_timeout && ~ask('continue?',options.ask_timeout)); b_stop = 1; end      
    
    if b_stop, 
        disp('stopping.'); 
        break; 
    end    
    nreassignments=[nreassignments(2:length(nreassignments)),ndif];
end;

% adjust pot
% % if done, the classes dont conform with pot!
% if ((ndif>0 && options.final_adjust) || options.remove_small && ~min(cl(:,end))) 
%     dprintf('calculating final potential and classes.');
%     nclend = max(cl(:,end));
%     suff = suffstat(nclend,type);
%     suff.adddata(data,cl(:,end));
%     pot = potential( suff, options.potential_options);
%     % TODO: finish
%     pr = freqs(cl(:,end))/sum(freqs(cl(:,end))); %TODO: or autropy
%     [ cl(:,end+1), potmin, potvals] = suff_data2v_TODO( data, pot, pr, options );
%     
%     freq_i = freqs(cl(:,end))';
%     freq(end+1,1:length(freq_i)) = freq_i;      
% end

dprintf('\ninitial & final class distributions:');
disp([freq(1,:);freq(end,:)]);

% assemble cl_info
if (max(cl) <= 255)
    cl_info.cl = uint8(cl);
else if max(cl) <= 255*255
    cl_info.cl = uint16(cl);        
    end
end

cl_info.freq = freq;
if (options.savedata)    
    cl_info.data = data;
end

cl_info.options = options;
cl_info.runtime = toc;
cl_info.potvals_priors = data2vpriors;
cl_info.potential_values = potvals;

dprintf('--------------end of suff_classify ------------\n\n');
end