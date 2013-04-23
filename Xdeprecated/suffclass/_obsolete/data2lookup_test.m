
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% data2lookup_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test data2lookup.m

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
     3     2     3
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
dfttrue=[0 0 1  4
         0 0 1 10
         1 2 5  2
         2 1 1  0];
[lookup,dft,conflict]=data2lookup(data,target);
conflict
input('next>');


load /pusers/neum/matlab/protein/DATA/bondRT.mat
ind=(max(alist,[],2)<=20);
data=alist(ind,:);
target=hlist(ind,:);
[lookup,dft,conflict]=data2lookup(data(1:N,:),target(1:N));
input('next>');

N=size(data,1);
[lookup,dft,conflict]=data2lookup(data(1:N,:),target(1:N));


