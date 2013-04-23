
*** not yet debugged ***
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
function [val,C]=value(s,i,acc,alp,lam);

count=shiftdim(s.count(i,1:s.type(i),:),1)'; % rather use permute!
nsum=sum(count);
mu=acc*count;
sig=sqrt(acc.^2*(count.*(1-count./nsum(ones(,:))));
mu=mu-alp*sig;

val=sum(max(mu))/sum(nsum);
C=mu'*mu+lam*(sig'*sig);





