% SUFF_CLASSIFY_TEST Test suff_classify and compare results with old code.
%
data=[
     1     1     3     4
     1     2     4     1
     1     3     3     4
     1     3     3     4
     1     4     2     1
     2     1     1     3
     2     1     3     1
     2     1     4     2
     2     2     1     3
     2     2     2     4
     2     2     3     1
     2     3     1     3
     2     3     4     2
     2     4     4     1
     3     1     1     4
     3     2     1     4
     3     2     3     2
     3     3     1     4
     3     3     2     1
     3     3     4     3
     3     4     3     1
     3     4     4     1
     4     2     1     1
     4     2     2     2
     4     2     3     3
     4     3     1     1
     4     3     1     1
     4     3     3     3
     4     4     1     1
     4     4     3     1];

target=data(:,4); % t=1 if y2=4, t=y1+y3 mod 4 otherwise 
data=data(:,1:3);

options.new_classes = 0;
clinfo=suff_classify(max(data,[],1),data,target,options);

addpath('_obsolete');
[lookup,pot,confusion]=catclass(data,target);
rmpath('_obsolete');

clsuff = clinfo.cl(:,end)'

clcc = lookup.g(data*lookup.fac - lookup.shift)'