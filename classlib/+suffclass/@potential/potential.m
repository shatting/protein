classdef potential < handle
% suffclass.potential < handle
%   Potential parameter estimation.
%
% PUBLIC PROPERTIES
%
%   .suff     suffclass.suffstat object used for creation
%   .options  options struct
%             .debug_usematcov  (0) use matlab to calculate numerical
%                                   moment
%             .num_reg_weight         weight of total potential
%                                               in regularization of
%                                               numerical class potentials
%                                               (1)
%             .cat_prior_reg          number of dummy categorical             
%                                               class points (0.01)
%             .cat_cond_reg           number of dummy categorical
%                                               conditional points
%                                               (sqrt(eps))
%             .cat_no_lookup          use cat lookup? (0)
%             .cat_lookup_saveclasspots     (0)
%
%   .numpot   numerical potential struct. see .calcnumpot()
%   .catpot   categorical potential struct. see .calccatpot()
%   .lookup   categorical lookup table
%   .matnumpot numerical potential struct using matlab covariance estimator
%
% CONSTRUCTORS
%
%   obj = potential(suffclass.suffstat suff)
%       Creates a suffclass.potential instance
%   obj = potential(suffclass.suffstat suff, additionaloptions)
%       Set additional options
%   obj = potential(suffclass.suffstat suff, additionaloptions || [], data, cl)
%       Use matlab numerical moment (options.debug_usematcov==1)
%       
% PUBLIC METHODS    
%
% .calcpot()
%   calculates .numpot and .catpot fields
%
% numpot = .calcnumpot( num_reg_weight ) (was SUFF_NUM2POT)
%   Generate numerical potential information. Sets .numpot.
%   INPUT: 
%       num_reg_weight(optional) 
%           regularization weight. 0 for no regularization. If >
%           0, total (over all classes) covariance estimation is added to 
%           each class covariance estimation with a weight of reg.
%
%   OUTPUT: (also sets .numpot property)
%       pot(cl)  potential information struct for class cl with fields:
%           .freq   sum of weights.
%           .mean   = sum(x)/.freq, mean estimation for this class.
%           .cov    = sum(x*x')/.freq - mean*mean', covariance estimation for this class.
%           .L      = (chol(.cov)')^-1, inverse lower cholesky factor of covariance matrix.
%           .ent    = log(det(.cov)), class entropy.
%
% catpot = .calccatpot( cat_prior_reg, cat_cond_reg ) (was SUFF_CAT2POT)
%   Generate categorical potential information from a sufficient statistic.
%   Sets .catpot.
%
%   INPUT:
%       cat_prior_reg     number of prior regularization points 
%       cat_cond_reg      number of conditionals regularization points
%
%   OUTPUT: (also sets .catpot property)
%       catpot  struct with fields:
%               .pot(i,a,cl)        potential for a in position i in class cl
%                   = -log(count(i,a,cl)+cat_cond_reg)+log(freq(c)+cat_cond_reg*length(.suff.cat));
%               .ent(cl)            cl-th class entropy    
%                   = -log(freq(cl) + cat_prior_reg);
%               .entsum             entsum(:,:,cl) == ent(cl). for simple
%                                   summation with .pot.
%                                   exp(-(pot.pot+pot.entsum)) ~ suff.count (~
%                                   b/c of regularization)
%               .freq               1 x ncl class freqs
%               .prob(i,a,cl)       P(C=cl|A_i = a)
%               .condprob(cl,a_1,..,a_d)     naive conditional probability
%                                   P( cl | A = (a_1,..,a_d) )
%       @smoothing: http://www.cs.cmu.edu/~tom/mlbook/NBayesLogReg.pdf    
%
% bool = .islookupdoable()
%   determines whether it is possible to calculate a categorical lookup
%   table (only possible if cat feature space is not too big)
%
% *lookup = calclookup( obj, potval_options ) (was SUFF_CATPOT2LOOKUP) (nyi)
%   Compute lookup table for a categorical potential. 
%   Done by iteration over all elements in feature space (slow).
%
%   lookup.v(y1,..,yd,c) := sum_i(pot.catpot(i,yi,c)) + entfactor*pot.ent(c)
%   lookup.g(y1,..,yd)   := argmin_c [ sum_i(pot.catpot(i,yi,c)) + entfactor*pot.ent(c) ]
%
%   Usage of lookup:
%   clpred = lookup.cl(data*lookup.fac-lookup.shift);
%       where data(l,:) is l-th feature vector
%       and clpred(l) is predicted class of l-th data vector
%
%   INPUT:
%       ymax(i)         number of cats in position i      
%       pot.pot(i,a,cl) potential for y_i=a in class cl
%       pot.ent(cl)     entropy of class cl
%       saveclasspots   boolean, all class potential values saved if true (default 0)
%       entfactor       (default 1)
%
%   OUTPUT:
%       lookup  .cl(y1,..,yd)       class of y
%               .v(y1,..,yd,cl)     potential value of y in class cl
%               .fac  
%               .shift
%
        
    properties
       suff
       options
       numpot
       catpot
       matnumpot       
    end
    
    methods           
        % constructor
        function obj = potential(suff, options, data, cl)
            obj.suff = suff;            
            if (nargin==1), options = struct; end
            obj.options = suffclass.potential.defaultoptions(options);
            
            obj.calcpot;
            
            if nargin == 4 && options.debug_usematcov,
                for cli = 1:max(cl),
                    cldata = data(cl==cli,suff.num);
                    obj.matnumpot(cli).cov = cov(cldata);
                    obj.matnumpot(cli).mean = mean(cldata)';
                    try
                        obj.matnumpot(cli).L = inv(chol(obj.matnumpot(cli).cov)');                        
                    catch exc
                        obj.matnumpot(cli).cov = obj.matnumpot(cli).cov.*(ones(length(suff.num))+diag(ones(1,length(suff.num)))*sqrt(eps));
                        obj.matnumpot(cli).L = inv(chol(obj.matnumpot(cli).cov)');                                                
                    end
                    obj.matnumpot(cli).ent = log(det(obj.matnumpot(cli).cov));
                    obj.matnumpot(cli).freq = size(cldata,1);
                end
                obj.numpot = obj.matnumpot;                
            end
        end             
                               
        function bool = islookupdoable(obj)
            bool = obj.suff.nnum == 0;
            try
                % see if lookuptable is generateable
                x = zeros(double(obj.suff.typ(obj.suff.cat)));
            catch
                bool = false;
            end
        end
        
    end
    
    methods (Access=protected)
        
        function calcpot( obj )
            
            % numeric portion
            if (obj.suff.nnum > 0),
                obj.calcnumpot(obj.options.num_reg_weight);
            end

            % categorical portion
            if (obj.suff.ncat > 0),    
                obj.calccatpot(obj.options.cat_prior_reg, obj.options.cat_cond_reg);
            end
        
        end
        
        function lookup = calclookup( obj, potval_options )
                              
            ncl = size(obj.catpot,3);            
            ymax = double(obj.suff.typ(obj.suff.cat));
            % initialize lookups; can throw error if too big
            try
                lookup.cl = zeros(ymax);                
            catch
                warning('not using cat lookup, feature space is too big (%i items)',prod(ymax));            
                lookup = [];
                return;
            end
            
            if obj.options.cat_lookup_saveclasspots,
                try    
                    lookup.v = zeros([ymax,ncl]); % could throw error even if .cl doesnt                    
                catch
                    lookup.v = [];
                    warning('not saving cat classpots.');                     
                end
            end

            % get fac and shift
            [~,lookup.fac,lookup.shift] = sys10toX(1,ymax);

            % run through all possible feature combinations
            % TODO: could be done faster&more elegantly i think
            if (0)
                casmax = prod(ymax(1:end));
                rep = fix(casmax/10);
                for cas = 1:casmax,
                   if (rem(cas,rep) == 0), dprintf('lookup %.2f%%',cas/casmax*100); end   

                   % get element
                   el = sys10toX(cas,ymax);

                   % get class potentials
                   v = suff_cat_vecpot(el, pot, entfactor);

                   % save class potentials for element
                    if saveclasspots
                       els = repmat(el,ncl,1);
                       iv = sysXto10([els [1:ncl]'],[ymax ncl]);
                       lookup.v(iv) = v;
                    end

                   % get minimizing class
                   [vs,lookup.cl(cas)] = min(v);
                end
            else
                casmax = prod(ymax(1:end));
                allcas = sys10toX(1:casmax,ymax);
                potvals = suffclass.potential_values(obj,potval_options);
                %TODO: finish
                %%potvals.calcvals(
            end
        end

        function numpot = calcnumpot( obj, num_reg_weight )

            moment = obj.suff.moment;
            empty = [];

            dim = size(moment,1) - 1;
            ncl = size(moment,3);

            if num_reg_weight > 0,
                momtot = sum(moment,3); 
                wtot = momtot(end,end);                
            end

            for c = 1:ncl,
                % regularization
                if num_reg_weight > 0,
                    moment(:,:,c) = moment(:,:,c) + (num_reg_weight/wtot)*momtot;    
                end
                numpot(c).freq = moment(end,end,c);

                if (numpot(c).freq == 0)        
                    numpot(c).mean = zeros(dim,1);
                    numpot(c).cov = eye(dim)*sqrt(eps);
                    numpot(c).L = eye(dim);
                    numpot(c).ent = inf;
                    empty = [empty, c];
                    continue;
                end
                numpot(c).mean = moment(end,1:end-1,c)' / numpot(c).freq; % nx1
                numpot(c).cov = moment(1:end-1,1:end-1,c);

                % form covariance
                numpot(c).cov = numpot(c).cov / numpot(c).freq - numpot(c).mean*numpot(c).mean';

                % regularize
                for i=1:dim,
                    numpot(c).cov(i,i)=numpot(c).cov(i,i)*(1+sqrt(eps));
                    if numpot(c).cov(i,i)==0, numpot(c).cov(i,i)=sqrt(eps); end;
                end

                try
                    numpot(c).L = inv(chol(numpot(c).cov)');
                catch exc
                    cov=numpot(c).cov
                    freq = numpot(c).freq
                    error('Class %i numerical covariance is not positive definite. Please use regularization.',c);                    
                end
                numpot(c).ent = log(det(numpot(c).cov));
            end

            if ~isempty(empty)
                warning('suffclass.potential.calcnumpot(): there were empty classes.');
                empty
            end
            
            obj.numpot = numpot;
        end
        
        function catpot = calccatpot( obj, cat_prior_reg, cat_cond_reg )

            if (nargin < 2)
                cat_prior_reg = 0;
            end
            if (nargin < 3)
                cat_cond_reg = 0;
            end
            dprintf('regularization with %g class points and %g conditionals.', cat_prior_reg, cat_cond_reg);

            count = obj.suff.count;
            ncl = obj.suff.ncl;
            ymax = obj.suff.typ(obj.suff.cat);
            
            d = length(ymax);
            maxcat = max(ymax);
            countdimmask = zeros(d,maxcat,ncl);
            for i=1:d, countdimmask(i,1:ymax(i),:) = 1; end;

            % counttot(i,a,c) = #category a is on place i, independent of c ~ P(a_i=a)
            counttot = sum(count,3);
            counttot = counttot(:,:,ones(ncl,1)); 

            % countclass(i,a,c) = #datapoints in class c
            countclass = repmat(sum(sum(count,1),2)/d,size(count,1),size(count,2));

            % .freq, class frequencies
            catpot.freq = squeeze(sum(count,2));
            catpot.freq = catpot.freq(1,:);

            % .pot
            if (cat_cond_reg == 0)
                ind = count>0;
                z = count==0;
                nzeroconds = sum(z(:));
                if (nzeroconds > 0)
                    dprintf('warning: number of zero conditional probabilities: %i. should use conditionals regularization.', nzeroconds);
                end
                catpot.pot = zeros(size(count)) + inf;
                catpot.pot(ind) = -log(count(ind)) + log(countclass(ind));
            else	
                dimensioncats = double(repmat(ymax',1,[maxcat ncl]));% strange error here with matlab2009b: second dimension is 127 long instead of 400
                %reason: ncl is int8, and matlab strangely converts the array to int8
                %and not to double (as it should i think)
                %fix: make ncl double in suffclass. no, make classes double    
                catpot.pot = -log(count+cat_cond_reg) + log(countclass.*dimensioncats+cat_cond_reg*dimensioncats);
                catpot.pot(~countdimmask) = inf;
            end

            % .ent
            if (cat_prior_reg == 0) 
               catpot.ent = zeros(1,ncl) + inf;
               indf = catpot.freq>0;
               nzeropriors = sum(catpot.freq==0);    
               if (nzeropriors > 0)
                   dprintf('warning: number of empty classes: %i. should use priors regularization.', nzeropriors);
               end
               catpot.ent(indf) = -log(catpot.freq(indf));
            else
               catpot.ent = -log(catpot.freq + cat_prior_reg);
            end

            % .entsum
            catpot.entsum = repmat(catpot.ent',1,[d maxcat]);
            catpot.entsum = permute(catpot.entsum,[2 3 1]);

            % .prob
            ind = counttot > 0;
            catpot.prob = zeros(d,max(ymax),ncl);
            catpot.prob(ind) = count(ind)./counttot(ind);
        
            obj.catpot = catpot;
        end            
        
    end
    
    methods (Static)
        function options = defaultoptions(options)
            
            if (nargin==0 || isempty(options))
                options = struct;
            end
                                 
            if ~isfield(options,'debug_usematcov'), options.debug_usematcov = 0; end;
                        
            if ~isfield(options,'num_reg_weight'),  options.num_reg_weight = 1; end;                        
            if ~isfield(options,'cat_prior_reg'),   options.cat_prior_reg = 0.01; end;
            if ~isfield(options,'cat_cond_reg'),    options.cat_cond_reg = sqrt(eps); end;

%             if ~isfield(options,'num_ent_weight'),  options.num_ent_weight = 1; end;
%             if ~isfield(options,'cat_ent_weight'),  options.cat_ent_weight = 1; end;
            
            if ~isfield(options,'cat_no_lookup'),   options.cat_no_lookup = 1; end;
            if ~isfield(options,'cat_lookup_saveclasspots'),   options.cat_lookup_saveclasspots = 0; end;
        end
    end

end

