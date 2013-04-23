
Contents file for 

matlab/ai/suffclass/
multivariate data visualization, clustering, and classification


ellipse.m        draw an ellipse
getfreq.m        calculate list of frequencies
ribbon.m         ribbon drawing of a chain of points
scatterplot.m    scatterplot of multivariate classified data
specclass.m      spectral line view of classified data vectors 
speclines.m      spectral lines representing a 1D vector 
specview.m       spectral line view of doubly stratified 1D data 


cat2freq.m
cat2group.m
catclass.m       partially supervised classification of categorical data
catcluster.m
catlookup.m
data2lookup.m

covclass.m       partially supervised classification of continuous data

cluster.m
clusterh.m       hierarchical clustering with Gini's measure
clusteri.m       iterative clustering with Gini's measure
clusters.m       sequential clustering with Gini's measure
dendro.m         dendrogram for hierarchical clustering
recluster.m      reclustering with Gini's measure

linkage.m	 hierarchical classification (by W. Huyer)


classdec.m       create decision tree 
                 *** pruning routine still missing!
classpot.m       create potential from sufficient statistic
                 classify data using the potential
mom2pot.m        create potential coefficients from moment information
pot2cl.m         classify data with a given potential
                 *** but classpot is more general!
suffadd.m        add sufficient statistics (for class*)
suffstat.m       handle sufficient statistic (for class*)
                 *** only opt=4 is implemented sofar!


value.m          assess classification quality of a categorical feature
                 *** not yet debugged


test programs
-------------

cat2freq_test.m
catcluster_test.m
catlookup_test.m
classdec_test.m
classpot_test.m
clusterh_test.m
clusters_test.m
cluster_test.m
data2lookup_test.m
dendro_test.m
mom2pot_test.m
recluster_test.m
scatterplot_test.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cat2freq.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function pairfreq=cat2freq(dft,lookup);
% function pairfreq=cat2freq(dft,lookup,0);
% create a pair frequency table pairfreq for categorical data 
% defined by data frequency table dft and dictionary lookup
%
% function pairfreq=cat2freq(data,target,1);
% create a pair frequency table pairfreq for categorical data 
% defined by labelled (data,target) pairs
%
% dft(y1,y2,...)        number of occurences of item (y1,y2,...)
% lookup                group dictionary; see cat2group.m
% data(l,:)             l-th data vector
% target(l,1)           l-th target (= class number to be predicted)
%
% pairfreq(yi,yk,g,ik)  frequency of y with y_i=yi,y_k=yk in group g
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cat2group.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function gpred=cat2group(lookup,data);
% function [gpred,confusion]=cat2group(lookup,data,target);
% function [gpred,confusion,g2t,mc]=cat2group(lookup,data,target,cost);
% classification of categorical data by table lookup
%
% data(l,:)             l-th data vector
% lookup                group dictionary
% lookup.g(y1,y2,...)   group label of item (y1,y2,...)
% lookup.fac            index factors
% lookup.shift          index shift
%                       dft(y1,y2,...)=dft(y*lookup.fac-lookup.shift)
% target(l,1)           l-th target (= class number to be predicted)
% cost(s,t)             cost of deciding for target s given target t
%                       must be nonnegative with zero diagonal
%                       (default: 1-eye, also used if cost=1)
%
% gpred(l,1)            predicted group for l-th data vector
% confusion(t,g)        number of targets t in group g; 
% g2t(1,g)=t            least costly class assignment t to group g
% mc(1,g)               misclassification cost for this assignment
% 
% total misclassification cost:    mctot=sum(mc);
% number of misclassified items:   nmisc=sum(g2t(gpred)~=target);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% catclass.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [lookup,pot,confusion]=catclass(data,target,dft,cost,costok);
% classification of categorical data defined by 
% labelled list (data,target) and unlabelled data frequency table dft
%
% data(l,:)             l-th data vector
% target(l,1)           l-th target (= class number to be predicted)
% dft(y1,y2,...)        number of unlabelled items (y1,y2,...)
%                       (default: uniform, also used if dft=1)
% cost(s,t)             cost of deciding for target s given target t
%                       must be nonnegative with zero diagonal
%                       (default: 1-eye, also used if cost=1)
% costok                acceptable total misclassification cost on data
%
% pot(yi,yk,g,ik)       potential for y with y_i=yi,y_k=yk in group g 
%                       ik=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 etc.
%                       i =1 1 2 1 2 3 1 2 3  4  1  2  3  4  5 
%                       k =2 3 3 4 4 4 5 5 5  5  6  6  6  6  6 
% lookup                group dictionary; see cat2group.m
% pot(yi,yk,g,ik)       revised potential 
% confusion(t,g)        number of targets t in group g; predict by
%                       tpred(g)=argmin(cost*confusion(:,g))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% catcluster.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [lookup,pairfreq,pot]=catcluster(dft,pot);
% clustering of categorical data defined by the data frequency table dft
% using an initial potential pot (or the desired number of clusters)
% based on alternating cat2class.m and cat2freq.m
%
% dft(y1,y2,...)        number of occurences of item (y1,y2,...)
% pot(yi,yk,g,ik)       potential for y with y_i=yi,y_k=yk in group g 
%                       ik=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 etc.
%                       i =1 1 2 1 2 3 1 2 3  4  1  2  3  4  5 
%                       k =2 3 3 4 4 4 5 5 5  5  6  6  6  6  6 
% pot=ng                only specifies the desired number of clusters
%                       (maximum 255)
% lookup                group dictionary; see cat2group.m
% pairfreq(yi,yk,g,ik)  frequency of y with y_i=yi,y_k=yk in group g
% pot(yi,yk,g,ik)       revised potential 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% catlookup.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [lookup,potsum]=catlookup(ymax,pot,report);
% create a dictionary lookup for categorical data of shape ymax
% with groups defined by the pair potential pot
%
% the dictionary is used by calling cat2group.m
%
% ymax(k)               number of alternatives in position k
% pot(yi,yk,g,ik)       potential for y with y_i=yi,y_k=yk in group g 
%                       ik=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 etc.
%                       i =1 1 2 1 2 3 1 2 3  4  1  2  3  4  5 
%                       k =2 3 3 4 4 4 5 5 5  5  6  6  6  6  6 
% report                (default inf) interval for showtime
%
% lookup                group dictionary; see cat2group.m
% potsum(y1,y2,...)     potential of item (y1,y2,...) in its group
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% classdec.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dectree,pot,pmin,cpred]=classdec(typ,data,N,in,out,options);
%          create decision tree data of type typ
%          using data(1:N,in) to predict data(1:N,out) 
%          rounding data to integers is assumed to be be harmless!! 
%          The potential pot (see classpot.m for details) allows one
%          to translate class information into data information
% function cpred=classdec(dectree,data,N);
%          classify incomplete data(1:N,:) using the decision tree
%          [this should not be used without preliminary pruning
%           but preferably postprocess with a pure potential method]
% function classdec(dectree);
%          print information about decision tree
%
% typ(i)     type of feature i
%            c>0: categorical variable with c levels
%            c<0: ordinal variable (meaningless mean) with |c| levels
%            c=0: numerical variable (meaningful mean)
% data(l,:)  l-th feature vector 
% N          number of vectors to be used
% in         list of variables on which decisions are based
% out        list of variables to be predicted
%            (may overlap; in=out for unsupervised classification)
% options    options influencing the tree construction
%            options(1)=opt defines recipe in potential training 
%	       opt = 1: C_i = C for all classes i
%	       opt = 2: C_i = sigma(i)*I (I unity matrix)
%	       opt = 3: C_i = diag(sigma(i,:))
%	       opt = 4: general C_i
%            options(2)=clsiz: minimal class size
%            options(3)=prc: percentage of data declared outliers
%            options(4)=prt: printlevel for classdec
%                       1: list decisions and times 
%                       2: also list loop through classes 
%            options(5)=cbysort: handling of categorical questions
%                       1: sort by min potential
%                       0: use simple clustering
%   
% cpred(l,1) class assignment of l-th feature vector
% dectree    information defining the decision tree
%            ncl=dectree(1)         number of classes 
%            e(node) entry point for node (e=2 at root)
%            i=dectree(e); status variable for node with entry point e
%            if i=0: leaf
%              cpred=dectree(e+1);  class assignment
%            if i>0: split by numerical or ordinal variable i 
%              e1=dectree(e+2);     entry point of child1
%              e2=dectree(e+3);     entry point of child2
%              sep=dectree(e+4);    separator for assignment to children
%                                   data(i)<sep: child1, else child2
%            for i<0: split by categorical, not ordinal variable -i
%              e1=dectree(e+1);     entry point of child1
%              e2=dectree(e+2);     entry point of child2
%              ee=dectree(e+3);     auxiliary index 
%              sep=dectree(e+4:ee); separator set for assignment
%                                   data(i) in sep: child1, else child2
% pot        information defining the potential for classes
%            (see classpot.m for details)
% pmin       potential minima at nodes; for pruning
%        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% classpot.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pot,potmin,freq]=classpot(suff,opt);
%          create potential from sufficient statistic 
%          rounding data to integers is assumed to be be harmless!! 
% function [cpred,prob]=classpot(pot,data,N,pr);
%          classify data(1:N,:) using the potential pot and a prior pr  
%
% suff       sufficient statistic
%            .ncl      number of classes
%            .typ      specifies type of variables
%            .ncat     number of categorical variables
%            .nnum     number of numerical variables
%            .nord     number of ordinal variables
%            .num      list of numerical features
%            .cat      list of categorical features
%            .ord      sublist of ordinal features
%                      ord(cat) is the list of ordinal features
%            .moment   symmetric moment matrix for numerical features
%                      moment(:,:,cl)= sum [x 1]^T*[x 1] 
%                      over all x=data(l,num) with class(l)=cl
%                      moment(end,end,:) is the class size
%                                        even if only categorical data
%            .count    count for categorical features
%                      count(i,a,cl) = number of l with 
%                                      class(l)=cl, data(l,cat(i))=a
% opt        defines recipe for potential creation (default; 4)
%            opt = 1: C_i = C for all classes i
% 	     opt = 2: C_i = sigma(i)*I (I unity matrix)
%	     opt = 3: C_i = diag(sigma(i,:))
%	     opt = 4: general C_i
% pot        information defining the potential
%            V(data) = Vnum(data(num)) + sum_i ccat(data(cat(i)))
%            Vnum(x) = (x-mu)'*R*(x-mu) + shift
%            the estimated conditional probability is
%            prob(data|cl)=exp(-V(data)) 
%            the total distribution needs a class prior
%            .ncl      number of classes
%            .typ      specifies type of variables
%            .ncat     number of categorical variables
%            .nnum     number of numerical variables
%            .cat      list of categorical features
%            .num      list of numerical features
%            .ccat     coefficients for categorical potential
%                      ccat(i,a,cl) for class cl, data(cat(i))=a
%            .mu       center vectors for numerical potential 
%                      mu(:,cl) for class cl
%            .sig      standard deviations for numerical potential
%                      sig(:,cl) for class cl
%            .R        quadratic term for numerical potential
%                      R(:,:,cl) (upper triangular) for class cl
%            .beta     scale factors for numerical potential
%            .shift    shifts for numerical potential
%                      shift(1,cl) for class cl
% data(l,:)  l-th feature vector (a row)
% N          number of vectors to be classified
% pr(cl)     prior probability or frequency of class cl
%            (default: all-one)
% cpred(l,1) predicted class assignment of l-th feature vector
% prob(l,cl) probability that l-th feature vector is in class cl
% potmin(cl) potential minimum of class cl 
% freq(cl)   regularized relative frequency of class cl in training
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% classpot_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test classpot.m, suffstat.m
% categorical/ordinal part not yet tested!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cluster.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [g2cl,mu,cl2c]=cluster(weight,sumdata,setting);
% clustering using Gini's measure 
% = vector quantization with constant scalar covariance
% calls the other clustering routines
%
% nrand random sequential clusterings 
% are followed by a consensus clustering
% relevant dendrograms are plotted
%
% weight(g)             number of vectors (or sum of weights) in group g
% sumdata(:,g)          (weighted) sum of vectors in group g
%                       sumdata(g,:) is allowed if sumdata is not square
% setting               parameter setting
%   .nrand                 number of random clusterings
%   .ncluster              number of initial clusters wanted
%   .ncnew                 number of clusters wanted after reclustering
%   .maxit                 maximal number of refinement iterations
%   .draw                  index set for drawing scatterplots
%   .ndraw                 maximal number of points drawn per cluster
%   .markersize            markersize for scatterplot
%
% g2cl(1,g)             final cluster to which group g is assigned
% mu(:,g)               weighted mean of vectors in final cluster g   
% cl2c(cl,nc)           class of cluster cl for a total of nc clusters
%                       (if fewer clusters are desired)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clusterh.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [merge,cost,D]=clusterh(countx,sumx);
% hierarchical clustering using Gini's measure 
% = vector quantization with constant scalar covariance
%
% function [merge,cost,sep]=clusterh(countx,sumx,splitlim);
% split a set into two by hierarchical clustering
% 
% countx(1,g)         number of vectors (or sum of weights) in group g
% sumx(:,g)           (weighted) sum of vectors in group g
% splitlim            minimal cluster weight for split (default 0)
%
% merge(:,i)=[l;k]    merge in step i group l into group k<l 
% cost(1,i)           cost of merging in step i
% D(f,g)              cost of merging initial groups f and g
% sep                 index list for a separating cluster 
%                     (next to root if splitlim<=1) 
%                     weight and complementary weight >=splitlim
%                     if this is impossible, sep=[]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clusteri.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [clabel,munew,cweight]=clusteri(weight,mudata,mu,maxit);
% iterative clustering using Gini's measure 
% = vector quantization with constant scalar covariance
%
% weight(1,g)         number of vectors (or sum of weights) in group g
% mudata(:,g)         (weighted) mean of vectors in group g
% mu(:,c)             center of cluster c
% maxit               maximal number of iterations (default; inf)
%
% clabel(1,g)         cluster label for group g
% munew(:,c)          center of groups labeled by cluster c
% cweight(1,c)        weight of cluster c
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clusters.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [weight,sumdata,D,mu]=clusters(weight,sumdata,ncl);
% sequential clustering using Gini's measure 
% = vector quantization with constant scalar covariance
%
% weight(1,g)         number of vectors (or sum of weights) in group g
% sumdata(:,g)        (weighted) sum of vectors in group g
% ncl                 number of final clusters 
% D(f,g)              cost of merging final clusters f and g
% mu(:,g)             weighted mean of vectors in final cluster g   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% covclass.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pot,cl,confus]=covclass(data,target,weight);
% partially supervised classification of continuous data 
% target labels are known classes if positive, tentative if negative, 
% unknown if 0
%
% The predicted classes may have different labels than initially.
% Predict original target labels by label=argmin(cost*confus(:,cl))
% where cost(s,t) is the cost of deciding for label s given target t
%
% data(l,:)             l-th data vector
% target(l,1)           l-th class number 
%                       fixed if positive, tentative if negative, 
%                       unknown if 0
% weight(l,1)           weight of data vector
% 
% pot                   potential defining classes by minimizing
%                       V_c(x) = ||L(x-mean)||^2+ent 
%                       over the classes with  V_c(x)<=Vmax_c 
%   .freq               class frequency
%   .mean               class mean 
%   .cov                class covariance
%   .ent                class entropy
%   .L                  class transform
%   .Vmax               maximal reliable potential value
% cl(l,1)               final class assignment
% confus(t,c)           number of labels +t in class c
%                       (but if all t<=0, of labels -t)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% data2lookup.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [lookup,dft,conflict]=data2lookup(data,target,opt);
% create a dictionary lookup and a data frequency table
% for categorical data given in (data,target) form
%
% data(l,:)             l-th data vector
% target(l,1)           l-th target (single class only)
% opt                   1: count conflicting data in dft (default)
%                       0: ignore conflicting data in dft
%                          (count only agreements with first)
%
% dft(y1,y2,...)        number of occurences of item (y1,y2,...)
% lookup                group dictionary; see cat2group.m
% conflict              conflict list: for l=conflict.
%                       data(l,:) has target distinct
%                       from corresponding lookup entry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dendro.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function g2c=dendro(weight,merge,cost,D,name);
% create dendrogramm (in the unit square) 
% for a merging scheme created by hierarchical clustering 
% (e.g., by clusterh.m)
%
% weight(1,g)         number of vectors (or sum of weights) in group g
% merge(:,i)=[l;k]    merge in step i group l into group k<l 
% cost(1,i)           cost of merging in step i
% D(f,g)              cost of merging groups f and g
%                     (if no preference, D=1);
% name(g,:)           name of group g (default: group number)
%
% g2c(g,nc)           class of group g for a total of nc clusters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ellipse.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h=ellipse(mu,r,C,lsty);
% draws an ellipse with center mu and covariance matrix r^2*C
% 
% lsty      line style (default '-b')
% h         line handle
%           set(h,'linestyle',sty)    changes the line style
%           set(h,'linewidth',wid)    changes the line width
%           set(h,'visible','off')    makes line invisible
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% freq2pot.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function pot=freq2pot(pairfreq);
% create potential from a pair frequency table
%
% pairfreq(yi,yk,g,ik)  frequency of y with y_i=yi,y_k=yk in group g
%
% pot(yi,yk,g,ik)       potential for y with y_i=yi,y_k=yk in group g 
%                       ik=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 etc.
%                       i =1 1 2 1 2 3 1 2 3  4  1  2  3  4  5 
%                       k =2 3 3 4 4 4 5 5 5  5  6  6  6  6  6 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% getfreq.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function freq=getfreq(list);
% calculate frequency table from a list of examples
% freq(a,b,...) is the number of i with list(i,1)=a, list(i,2)=b, ...
%
% function freq=getfreq(list,wt);
% calculate frequency table from a weighted list of examples
% freq(a,b,...) is the sum of wt(i) with list(i,1)=a, list(i,2)=b, ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% linkage.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ord,c,bounds,diss] = linkage(d,opt)
% hierarchical cluster analysis by single or complete linkage 
% 
% Input:
% d		symmetric matrix of dissimilarities between classes i 
%		and j
% opt 		opt=0: single linkage, opt=1: complete linkage
%
% Output:
% ord		permutation of the set of classes such that the fusions 
%		in the clustering algorithm can be performed without 
%		crossing lines
% c		In the k-th step of the cluster analysis, two clusters
% 		with dissimilarity c(k) are joined
% bounds	ord(bounds(1,k):bounds(2,k)) gives the elements of the
%		cluster generated in the k-th step
% diss		vector of length n-1 of dissimilarities between 
%		joining clusters; joining classes ord(j), ord(j+1) 
%		when diss(j) <= c(i) gives n - i clusters
%		c = sort(diss)
%		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mom2pot.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [meanx,Rx,beta,Vmin,sigx]=mom2pot(countx,sumx,Mx,reg);
% function [meanx,Rx,beta,Vmin,sigx]=mom2pot(suff,reg);
%
% create potential coefficients from moment information
% generated elsewhere from data vectors with a class assignment
%
%
% function mom2pot(prt);  % for debug
% repeat calculation with previous data and printlevel prt
% 
% Important: means are rounded to integers; scale data if necessary! 
% 
% The potential has the form
%    V(x,t)=beta(t)*(x-meanx(:,t))'*Rx(:,:,t)*(x-meanx(:,t))+Vmin(t)
% with integral Rx, and is related to the Gaussian densities by
%    p(x|t)=exp(-V(x,t)).
% For prediction with prior class probabilities pr(t),
% minimize V(t|x)=V(x,t)-log(pr(t))
%
% The covariance matrices can be obtained by
%    C=inv(Rx(:,:,t)+Rx(:,:,t)')/beta(t);
% the potential mean is
%    Vmean(t)=Vmin(t)+nd/2
% where nd=size(sumx,1) is the dimension of the data vectors.
%
% Regularization ensures that R+R' is positive definite 
% 
% countx(1,t)     number of vectors of class t
% sumx(:,t)       sum  of vectors of class t
% Mx(:,:,t)       2nd moment matrix of vectors of class t
% suff            sufficient statistic containing countx,sumx,Mx
% reg             number of dummy entries, must be >1 (default 2)
%
% meanx(:,t)      mean for class t (rounded to integer)
%                 or empty (if data set is too small; 
%                           then Rx contains the error message)
% sigx(:,t)       standard deviations for class t (optional)
% Rx(:,:,t)       upper triangular quadratic term for class t
% beta(1,t)       potential scaling factor
% Vmin(1,t)       potential minimum
%                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pot2cl.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pred,ppot]=pot2cl(meanx,Rx,beta,Vmin,data,logpr);
% predict classes for data from potential
%    V(x,t)=beta(t)*(x-meanx(:,t))'*Rx(:,:,t)*(x-meanx(:,t))+Vmin(t)
% created by mom2pot.m
%
% For prediction with prior class probabilities pr(t),
% minimize V(t|x)=V(x,t)-log(pr(t)),
% i.e., use Vmin-log(pr) in place of Vmin
%
% meanx(:,t)      mean for class t (rounded to integer)
% Rx(:,:,t)       upper triangular quadratic term for class t
% beta(1,t)       potential scaling factor
% Vmin(1,t)       potential minimum
% data(l,:)       l-th data vector, = x^T
%
% pred(l,1)       predicted class of l-th data vector
% ppot(l,1)       selected potential of l-th data vector
%                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% recluster.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [g2cl,mu,cl2c]=recluster(weight,sumdata,g2cl,maxit,ncnew);
% reclustering using Gini's measure 
% = vector quantization with constant scalar covariance
% calls the other clustering routines
%
% weight(g)             sum of weights in group g
% sumdata(:,g)          weighted sum of vectors in group g
% g2cl(1,g)             initial cluster to which group g is assigned
%                       (if too many clusters, use insted cluster.m)
% maxit                 maximal number of refinement iterations
% ncnew                 merge to ncnew clusters (or fewer)
%
% g2cl(1,g)             final cluster to which group g is assigned
% mu(:,g)               weighted mean of vectors in final cluster g 
% cl2c(cl,nc)           class of cluster cl for a total of nc clusters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ribbon.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ribbon(x,scl,thr);
% makes six ribbon drawings of a 3D curve obtained by linear 
% interpolation of x(1,:), ..., x(n,:)
% scl:  scale factor for ribbon width (default 1)
% thr:  number of threads (default 5)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% scatterplot.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function scatterplot(data,class,classname,compname,ndraw, markersize);
% draw a thinned-out scatterplot with up to ndraw points per class
% for multidimensional data, all pair projections are drawn
%
% data(l,:)	  l-th data point
% class(l,1)	  class of l-th data point
%                 but class=1 for single class
% classname(c,:)  name of class c  (or '' for default c)
% compname(k,:)   name of component k (or '' for default k)
% ndraw           maximal number of points per class drawn
% markersize      marker size (1=smallest)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% specclass.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,y]=specclass(data,class,classname,compname)
% spectral line view of a collection of classified data vectors 
% (many rows are meaningful, but not many columns)
%
% data(l,:)	  l-th data point
% class(l)	  class of l-th data point
% classname(c,:)  name of class c
% compname(k,:)   name of component
%
% x,y		  vectors for reconstructing the plot via
%                    for i=1:size(x,1),
%                      plot(x(i,:),y(i,:),'k-');
%                    end;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% speclines.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function speclines(x,nr,range,spectext,sparse)
% prints spectral lines representing a 1D vector x
% with range range and ordinate [nr,nr+0.8] 
% text in spectext (default: none) is appended after the spectrum
% if sparse (default sparse=0), a random point on each spectral line 
%    is drawn instead of the whole spectral line, to see the density 
%    when there are too many points 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% specview.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [sizes,x,y]=specview(data,row,col,rowname,colname)
% spectral line view of a collection of 1D data 
% stratified by two categorical variables
% (many rows are meaningful, but not many columns)
%
% data(l)	l-th data point
% row(l)	row containing l-th data point
% col(l)	column containing l-th data point
%		(or 1 for single column)
% rowname(r,:)  name of row r
% colname(c,:)  name of column c
%
% sizes(r,c)	number of data points in cell (r,c)
% x,y		vectors for reconstructing the plot via
%                  for i=1:size(x,1),
%                    plot(x(i,:),y(i,:),'k-');
%                  end;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% suffadd.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function suff=suffadd(suff1,suff2,cas);
% add sufficient statistics suff1, suff2 and assign it to suff
% cas=1:     same classes
% cas=2:     create new class for each suff2 class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% suffstat.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function suff=suffstat(ncl,typ);
%          initialize sufficient statistic 
%          for data of arbitrary type (>1 feature) and ncl classes
% function suff=suffstat(suff,data,N,clist);
%          increment sufficient statistic using data(1:N,:) 
%          with assigned classes clist(1:N)
% function suff=suffstat(suff,data,N,clist,cold);
%          modify sufficient statistic using data(1:N,:) when
%          class assignments change from cold(1:N) to clist(1:N)
% function suff=suffstat(suff,ncl);
%          increase the number of classes to ncl
%
% ncl        number of classes
% typ(i)     type of feature i
%            c>0: categorical variable with c levels
%            c<0: ordinal variable (meaningless mean) with |c| levels
%            c=0: numerical variable (meaningful mean)
% suff       sufficient statistic
%            .ncl      number of classes
%            .typ      specifies type of variables
%            .ncat     number of categorical variables
%            .nnum     number of numerical variables
%            .nord     number of ordinal variables
%            .num      list of numerical features
%            .cat      list of categorical features
%            .ord      sublist of ordinal features
%                      ord(cat) is the list of ordinal features
%            .moment   symmetric moment matrix for numerical features
%                      moment(:,:,cl)= sum [x 1]^T*[x 1] 
%                      over all x=data(l,num) with class(l)=cl
%                      moment(end,end,:) is the class size
%                                        even if only categorical data
%            .count    count for categorical features
%                      count(i,a,cl) = number of l with 
%                                      class(l)=cl, data(l,cat(i))=a
% data(l,:)  l-th feature vector in current batch
% N          number of vectors in current batch
% clist(l,1) class assigned to l-th feature vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% unitest.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [fac,perc,kuiper]=unitest(data);
% Kuiper test statistics for testing uniformity of data 
%
% roughly perc percent of uniform samples exceed the test statistic
% fac   1.7  1.8  1.9  2.0  2.1  2.2   2.3   2.4   2.5   2.6   2.7   
% perc  6.8  4.6  3.1  1.8  1.1  0.56  0.28  0.11  0.06  0.03  0.015
%
% thus fac<=1.8, giving
%   kuiper<1.8*(sqrt(N-1.5)-0.7),
% is a good uniformity test
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% value.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [val,C]=value(suff,i,acc,alp,lam);
% assess the classification quality of a categorical feature
%
% suff       sufficient statistic
%            .ncl      number of classes
%            .ncat     number of categorical variables
%            .cat      list of categorical features
%            .count    count for categorical features
%                      count(i,a,cl) = number of l with 
%                                      class(l)=cl, data(l,cat(i))=a
% i         feature whose value is to be found
% acc       acceptability matrix
%           c1 is acceptable as a prediction for c2 
%           if acc(c1,c2) is close to 1
% alp       safety factor (default 1)
% lam       weight of std contribution (default 0 if alp>0, 1 if alp=0)
% val       value of feature i
%           features with values close to 1 are highly predictive
% C         moment matrix of feature i
%
