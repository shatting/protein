% DENDRO_TEST Test dendro.m



% function dendro(countx,merge,cost,D,name);
% create dendrogramm for a merging scheme 
% created by hierarchical clustering (e.g., by clusterh.m)
%
% countx(1,g)         number of vectors (or sum of weights) in group g
% merge(:,i)=[l;k]    merge in step i group l into group k<l 
% cost(1,i)           cost of merging in step i
% D(f,g)              cost of merging groups f and g
% name(g,:)           name of group g (default: group number)
%

% test example
ng=8;
countx=[7 8 20 5 4 6 26 24];
merge=[2 3 6 8 7 4 5
       1 1 5 7 5 1 1];
cost= [1 2 3 4 5 6 7];
D=ones(8);
D=rand(8);
% load D
name=[]; % use group number
g2c=dendro(countx,merge,cost,D,name)
