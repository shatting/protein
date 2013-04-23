

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% cat2freq_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test cat2freq.m 

data=[
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
     4     2     2
     3     2     2
     4     1     2
     4     1     2
     3     1     2
     3     3     3
     3     3     3
     3     3     3
     3     4     3
     4     3     3
     3     4     3
     3     3     3
     2     3     3
     3     2     3
     3     3     3
     1     4     4
     1     4     4
     1     3     4
     1     4     4
     1     4     4];
data=[data,data];
d=size(data,2);
target=sum(data,2)-5; % uniquely determined by data
ng=max(target);
if min(target)<1,
  mintarget=min(target),
  error('bad target'); 
end;

[lookup,dft,conflict]=data2lookup(data,target);
conflict % should be empty

% pairfreq(yi,yk,g,ik)  frequency of y with y_i=yi,y_k=yk in group g
pairfreq=cat2freq(dft,lookup);
pairfreq1=cat2freq(data,target,1);
for g=1:ng,
  for ik=1:6,
    if sum(sum(pairfreq(:,:,g,ik)~=(pairfreq1(:,:,g,ik))))>0,
      g,ik
      freq=[pairfreq(:,:,g,ik),pairfreq1(:,:,g,ik)]
      error('programming error');
    end;
  end;
end;

disp('*** if this works, repeat with data=[data,data] ***');




