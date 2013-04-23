%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% cat2group_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test cat2group.m

% data(l,:)             l-th data vector
% lookup                group dictionary; 
% lookup.g(y1,y2,...)   group label of item (y1,y2,...)
% lookup.fac            index factors
% lookup.shift          index shift
%                       dft(y1,y2,...)=dft(y*lookup.fac-lookup.shift)
% target(l,1)           l-th target (= class number to be predicted)
% cost(s,t)             cost of deciding for target s given target t
%                       (default: 1-eye, also used if cost=1)
%
% gpred(l,1)            predicted group for l-th data vector
% confusion(t,g)        number of targets t in group g; 
% g2t(1,g)=t            least costly class assignment t to group g
% mc(1,g)               misclassification cost for this assignment
% 
% total misclassification cost:    mctot=sum(mc(gpred));
% number of misclassified items:   nmisc=sum(g2t(gpred)~=target);
%
[gpred,confusion,g2t,mc]=cat2group(lookup,data,target,cost);


A=[
     1     3     4
     1     4     4
     1     4     4
     1     4     4
     1     4     4
     2     3     3
     2     4     1
     2     4     1
     2     4     1
     2     4     1
     2     4     1
     2     4     1
     2     4     1
     2     4     1
     2     4     1
     2     4     1
     3     1     2
     3     2     2
     3     2     2
     3     3     3
     3     3     3
     3     3     3
     3     3     3
     3     3     3
     3     4     3
     3     4     3
     4     1     2
     4     1     2
     4     2     2
     4     3     3];
data=A(:,1:2);target=A(:,3);
[lookup,dft,conflict]=data2lookup(data,target);
% modify lookup
lookup.g=[
     0     0     4     4
     0     0     3     1
     2     2     3     3
     2     2     3     0]; % true lookup
lookup.g=[
     1     0     4     2
     0     0     3     1
     2     2     3     3
     2     2     4     0]; % modified lookup
[gpred,confusion,g2t,mc]=cat2group(lookup,data,target,1);
confusion,mc
miscc=sum(gpred~=target)

