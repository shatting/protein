% SUFFCLASS
% A MATLAB package for data analysis.
%
% This package is intended for use in fully or partial classification and
% cluster analysis of data consisting of real ('numerical') and/or integer
% (categorical (meaningless mean) and ordinal (meaningful mean)) features. 
% 
% Usage: First, run startup.m. 
% The basic data structure is generated with the use of suffstat.m, an 
% implementation of the notion of sufficient statistics. Data then can be 
% added to the structure using the same function. Then, a potential 
% estimation ('potential') can be built using suff_pot. Notice however, 
% that such a 'potential' is in fact a collection of a number of potentials, 
% each stemming from one class defined in the data adding step. Simple 
% classification can then be done with the use of suff_vecpot (for single 
% vectors), suff_data2cl (which is essentially iterating suff_vecpot) or 
% the fast suff_data2v (using MATLAB vectorization).
% 
% There are, however, several algorithms for convenient and more elaborate 
% classification and clustering tasks:
%
% suff_classify: partially or fully supervised iterative classification.
%   This is basically iterating the steps described above (suff. statistic
%   ->potential estimation->reclassification). New classes can be
%   introduced automatically via one classification step's confusion
%   matrix (large entries off the diagonal get a new class assignment) and
%   small classes can be removed automatically.
%
% suff_clusterh: hierarchical cluster analysis. Given classes are
%   sucessively merged based on lowest cost in total potential value.
%   cluster_tryfindcut then can be used to automatically find a good cut in
%   the merge sequence obeying desired parameters (eg. number of classes).
% 
%
% Files
%   startup                    - startup script for suffclass package.
%   applymerge                 - Apply merging information to a list of class numbers.
%   cluster_plot               - Dendrogramm and merge cost plot, visualization of cut finding criteria.
%   cluster_tryfindcut         - Try to find a good cut in cluster merge sequence.
%   suff_cat2pot               - Generate categorical potential information from a sufficient statistic.
%   suff_cat_vecpot            - Calculate potential value of a categorical feature vector.
%   suff_catpot2lookup         - Compute lookup table for a categorical potential. 
%   suff_classify              - Partially or fully supervised iterative classification of mixed data.
%   suff_classify_test         - Test suff_classify and compare results with old code.
%   suff_cluster               - Cluster analysis, cut finding and applying.
%   suff_clusterh              - Hierarchical cluster analysis of (only categorical yet) data.
%   suff_data2cl               - Non-vectorized class prediction.
%   suff_data2v                - Vectorized class prediction and potential evaluation.
%   suff_num2pot               - Generate numerical potential information from a sufficient statistic.
%   suff_pot                   - Generate potential information from sufficient statistic.
%   suff_pottest               - Test suff_data2cl and suff_data2v and compare classifiaction results.
%   suff_vecpot                - Calculate potential value of a (mixed) feature vector.
%   suffadd                    - Add two sufficient statistics.
%   suffstat                   - Create, add and remove data to or from a sufficient statistic.
%   scr_compare                - Test and compare old and new potential estimators.
%
% *************************************************************************
% 
% APPLYMERGE Apply merging information to a list of class numbers.
%  
%   [ clmerged, m ] = applymerge( cl, merge )
%   INPUT:
%     cl(l,1)     l-th class
%     merge       in ith step, merge merge(1,i) into merge(2,i)
%  
%   OUTPUT: 
%     clmerged    = m(cl)
%     m           merge information after all steps
%  
%   See also: suff_clusterh.
% 
% CLUSTER_PLOT Dendrogramm and merge cost plot, visualization of cut finding criteria.
%   cluster_plot( cl, cluster_info, cluster_result, options )
% 
% CLUSTER_TRYFINDCUT Try to find a good cut in cluster merge sequence.
%   [ cut ] = cluster_tryfindcut( cl, cost, options)
% 
% SCR_COMPARE Test and compare old and new potential estimators.
%   Example usage: 
%     set new=1, alg=1 and run. set new=0, alg=2, and run again. 
% 
% STARTUP startup script for suffclass package.
% 
% SUFF_CAT2POT Generate categorical potential information from a sufficient statistic.
% 
%   catpot = SUFF_CAT2POT( suff, rpriors, rconditionals ) 
%  
%   catpot.pot(i,a,cl) = -log(count(i,a,cl)+rconditionals)+log(freq(c)+rconditionals*d);
%   catpot.ent(cl)     = -log(freq(cl) + rpriors);
%   (d = length(suff.cat))
%  
%   INPUT:
%     suff        sufficient statistic
%     rpriors     number of prior regularization points 
%     rcond       number of conditionals regularization points
%  
%   OUTPUT:
%     catpot      .pot(i,a,cl)        potential for a in position i in class cl
%                 .ent(cl)            cl-th class entropy    
%                 .entsum             entsum(:,:,cl) == ent(cl). for simple
%                                     summation with .pot.
%                                     exp(-(pot.pot+pot.entsum)) ~ suff.count (~
%                                     b/c of regularization)
%                 .freq               1xncl class freqs
%                 .prob(i,a,cl)       P(C=cl|A_i = a)
%  
%   @smoothing:
%   http://www.cs.cmu.edu/~tom/mlbook/NBayesLogReg.pdf    
%  
%   See also: suffstat.
% 
% SUFF_CAT_VECPOT Calculate potential value of a categorical feature vector.
% 
%   This function is very slow.
%  
%   v = suff_cat_vecpot( x, catpot, entfactor )
%  
%   INPUT: 
%     x           1xd categorical feature vector.
%     catpot      categorical potential (see suff_cat2pot)
%     entfactor   factor to apply to potential entropies.
%  
%   OUTPUT:
%     v           nclx1 class potential values.
%  
%   See also: suff_vecpot, suff_cat2pot.
% 
% SUFF_CATPOT2LOOKUP Compute lookup table for a categorical potential. 
% 
%   Done by iteration over all elements in feature space (slow).
%  
%   lookup = SUFF_CATPOT2LOOKUP( ymax, pot, saveclasspots, entfactor )
%  
%   lookup.v(y1,..,yd,c) := sum_i(pot.catpot(i,yi,c)) + entfactor*pot.ent(c)
%   lookup.g(y1,..,yd)   := argmin_c [ sum_i(pot.catpot(i,yi,c)) + entfactor*pot.ent(c) ]
%  
%   Usage of lookup:
%   clpred = lookup.cl(data*lookup.fac-lookup.shift);
%     where data(l,:) is l-th feature vector
%     and clpred(l) is predicted class of l-th data vector
%  
%   INPUT:
%     ymax(i)         number of cats in position i      
%     pot.pot(i,a,cl) potential for y_i=a in class cl
%     pot.ent(cl)     entropy of class cl
%     saveclasspots   boolean, all class potential values saved if true (default 0)
%     entfactor       (default 1)
%  
%   OUTPUT:
%     lookup  .cl(y1,..,yd)       class of y
%             .v(y1,..,yd,cl)     potential value of y in class cl
%             .fac  
%             .shift
%  
%   See also: suff_catpots.
% 
% SUFF_CLASSIFY Partially or fully supervised iterative classification of mixed data.
% 
%   cl_info = suff_classify(type, data, target, options)
%  
%   target labels are known classes if positive, unknown if 0
%  
%   The predicted classes may have different labels than initially.
%   Predict original target labels by label=argmin(cost*confus(:,cl))
%   where cost(s,t) is the cost of deciding for label s given target t
%  
%   INPUT:
%   type(i)               type of feature i (see suffstat.m)
%   data(l,i)             i-th feature of datapoint l
%   cl(l,1)               initial class assignment
%                         unknown if 0
%   options               .plot                   should plots be done (0)
%                         .max_iterations         maximum iterations (inf)
%                         .ask                    ask after each iteration(0)
%                         .ask_timeout            timeout for ask(5)
%                         .term_percentage        terminate if
%                                                 relabeled<this (0.05)
%                         .new_classes            what percentage of
%                                                 total data points should be
%                                                 exceeded in iteration
%                                                 confusion matrix to form
%                                                 new class (0.025)
%                         .remove_small           
%                         .stop_on_nonewclasses   stop if no new classes are
%                                                 formed (0)
%                         .num_reg_weight         weight of total potential
%                                                 in regularization of
%                                                 numerical class potentials
%                                                 (1)
%                         .num_ent_weight         weight of entropy term in 
%                                                 numerical part (1)
%                         .cat_prior_reg          number of dummy categorical             
%                                                 class points (0.01)
%                         .cat_cond_reg           number of dummy categorical
%                                                 conditional points
%                                                 (realmin)
%                         .cat_ent_weight         weight of categorical
%                                                 priors (0.5). 0 for uniform
%                         .cat_no_lookup          use cat lookup? (0)
%                         .final_adjust           adjust potential after
%                                                 stopping to reflect final
%                                                 class changes? (0)
%                         .savedata               data saved to output? (0)
%                                                 .
% 
% SUFF_CLASSIFY_TEST Test suff_classify and compare results with old code.
% 
% SUFF_CLUSTER Cluster analysis, cut finding and applying.
% 
% SUFF_CLUSTERH Hierarchical cluster analysis of (only categorical yet) data.
% 
%   [merge,cost,Dsuff0,Dsuffend]=suff_clusterh(type, data, cl, options, verbose)
%   hierarchical clustering of row vectors contained in data, having features described by type
%   and preassigned classes cl.
%  
%   [merge,cost,Dsuff0,Dsuffend]=suff_clusterh(suff, options, verbose)
%   hierarchical clustering of the classes described by suff
%   
%   INPUT:
%  
%   type(i)             type of feature i
%   data(l,:)           l-th data vector
%   cl(l)               class of l-th data vector
%   
%   suff                suffictient statistic describing the data and classes, see suffstat.m 
%  
%   options             not used as of now
%  
%   OUTPUT: 
%  
%   merge(:,i)=[l;k]    merge in step i class l into class k<l 
%   cost(1,i)           cost of merging in step i
%   Dsuff0(f,g)         initial merge cost matrix
%   Dsuffend(f,g)       final cost matrix
%  
%   See also: suff_classify.
% 
% SUFF_DATA2CL Non-vectorized class prediction.
% 
%   clpred = SUFF_DATA2CL( data, pot, options )
%  
%   Calculate predicted class assignments for data under potential pot. Only
%   works fast for categorical features with a lookup table.
%   Otherwise, iterates suff_vectorpotential.
%   For a vectorized version, see suff_data2v.
%  
%   INPUT:
%     data        nxd array, n feature vectors of dimension d
%     pot         potential struct (see suff_pot)
%     options     
%  
%   See also: suff_data2v, suff_pot, suff_vectorpotential.
% 
% SUFF_DATA2V Vectorized class prediction and potential calculation.
% 
%   [ cl, vmin, V, prob ] = SUFF_DATA2V( data, pot, pr, options )
%   Calculate class assignment, potential values and conditional
%   probabilities for data under potential pot. This is the vectorized
%   version of suff_data2cl.
%  
%   INPUT:
%   pot           potential (suff_pot.m)
%   data(l,:)     l-th data vector
%   pr            (optional) input priors. if not given, then uniform assumed
%   options       .ent_weight (default 1)
%                 .use_entropy (default false)    if true, pr is assumed to
%                                                 be entropy values instead
%                                                 of probabilities.
%  
%   OUTPUT: 
%   cl(l,1)       predicted class for lth data vetor
%   vmin(l,1)     minimum potential value
%   V(l,cl)       potential value of lth vector under potential cl
%   prob(l,cl)    probability that lth vector is in class cl
%  
%   See also: suff_pot, suff_data2cl.
% 
% SUFF_NUM2POT Generate numerical potential information from a sufficient statistic.
%
% pot = suff_num2pot( suff, reg )
% 
% INPUT: 
%   suff    Sufficient statistic.
%   reg     (optional) regularization weight. 0 for no regularization. If >
%           0, total (over all classes) covariance estimation is added to 
%           each class covariance estimation with a weight of reg.
%
% OUTPUT:
%   pot(cl)  potential information struct for class cl with fields:
%       .freq   number of points.
%   	.mean   = sum(x)/.freq, mean estimation for this class.
%   	.cov    = sum(x*x')/.freq - mean*mean', covariance estimation for this class.
%   	.L      = (chol(.cov)')^-1, inverse lower cholesky factor of covariance.
%   	.ent    = log(det(.L)), class entropy.
%
% See also: suff_cat2pot, suff_pot.
%
% SUFF_POT Generate potential information from sufficient statistic.
% 
%   pot = SUFF_POT ( suff, options )
%   create potential from sufficient statistic.
%  
%   INPUT:
%     suff        sufficient statistic (see suffstat)
%     options     (optional) 
%                 .cat_cond_reg = realmin
%                 .cat_prior_reg = 0.1
%                 .cat_ent_weight = 0.5
%                 .cat_no_lookup = 1
%                 .reg_weight = 1
%  
%   OUTPUT:
%    pot      potential struct with fields
%                 .numpot     see suff_num2pot.m    
%                 .catpot     see suff_cat2pot.m
%  
%   See also: suffstat, suff_num2pot, suff_cat2pot.
% 
% SUFF_POTTEST Test suff_data2cl and suff_data2v and compare classifiaction results.
% 
% SUFF_VECPOT Calculate potential value of a (mixed) feature vector.
% 
%   [ v, cl, V ] = SUFF_VECPOT( x, pot, options )
%   Calculate the (minimum) potential value(s) of a feature vector.
%  
%   INPUT:
%     x       the feature vector
%     pot     the potential (@see suff_pot.m)
%     options (optional)
%             .cat_ent_weight (default 1)
%             .num_ent_weight (default 1)
%   OUTPUT:
%     [v, cl] = min(V) and V(c) = pot.cat(x,c) + pot.num(x,c) for all c
%   
%   See also: suff_data2cl, suff_data2v.
% 
% SUFFADD Add two sufficient statistics.
% 
%   function suff=suffadd(suff1,suff2,cas);
%   add sufficient statistics suff1, suff2 and assign it to suff
%   cas=1:     same classes
%   cas=2:     create new class for each suff2 class
% 
% SUFFSTAT Create, add and remove data to or from a sufficient statistic.
% 
%   function suff=suffstat(ncl,typ);
%            initialize sufficient statistic 
%            for data of arbitrary type (>1 feature) and ncl classes
%  
%   function suff=suffstat(suff,data,N,clist,[],wt);
%            increment sufficient statistic using data(1:N,:) 
%            with assigned classes clist(1:N) and weights wt
%  
%   function suff=suffstat(suff,data,N,clist,cold,wt,wtold);
%            modify sufficient statistic using data(1:N,:) when
%            class assignments change from cold(1:N) to clist(1:N)
%            with new weights wt and old weights wtold
%  
%   function suff=suffstat(suff,ncl);
%            increase the number of classes to ncl
%  
%   function suff=suffstat(suff,mergelist);
%            merge classes from mergelist 
%  
%   ncl        number of classes
%   typ(i)     type of feature i
%              c>0: categorical variable with c levels
%              c<0: ordinal variable (meaningless mean) with |c| levels
%              c=0: numerical variable (meaningful mean)
%   suff       sufficient statistic
%              .ncl      number of classes
%              .typ      specifies type of variables
%              .ncat     number of categorical variables
%              .nnum     number of numerical variables
%              .nord     number of ordinal variables
%              .num      list of numerical features
%              .cat      list of categorical features
%              .ord      sublist of ordinal features
%                        ord(cat) is the list of ordinal features
%              .moment   symmetric moment matrix for numerical features
%                        moment(:,:,cl)= sum [x 1]^T*[x 1] 
%                        over all x=data(l,num) with class(l)=cl
%                        moment(end,end,:) is the class size
%                                          even if only categorical data
%              .count    count for categorical features
%                        count(i,a,cl) = number of l with 
%                                        class(l)=cl, data(l,cat(i))=a
% 
