
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% catclass_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test catclass.m

clear;
tic;
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
     3     3     1     1
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
d=size(data,2);
ng=max(target)
dft=1;   % uniform distribution
cost=-1; % unit gain for agreement
data=uint8(data);target=uint16(target);
disp('before catclass');
[lookup,pot,confusion]=catclass(data,target,dft,cost);

input('next>');
data=double(data);
data=[data,sum(fix((data*2/3)),2)]; % and old target
ng=max(target)
if min(target)<1,
  mintarget=min(target),
  error('bad target'); 
end;
dft=1;   % uniform distribution
cost=-1; % unit gain for agreement
disp('before catclass');
[lookup,pot,confusion]=catclass(data,target,dft,cost);

input('next>');

% partially supervised secondary structure prediction HHHHH vs. EEEEE
load bondRT.mat
ind=(max(alist,[],2)<=20 & (hlist==1|hlist==119));
alist=alist(ind,:);
hlist=hlist(ind);
hlist(hlist==1)=1;
hlist(hlist==119)=2;
hlist(hlist==60)=3;
cost=1; % unit gain for agreement

ind=1:500;                   % labelled part of the data
fragind=[2:4];               % part used as fragments
data=alist(:,fragind);
target=hlist;
[N,d]=size(data)
disp('before data2lookup');
[lookup,dft,conflict]=data2lookup(data,target);
nconflict=length(conflict)
disp('before catclass');
[lookup,pot,confusion]=catclass(data(ind,:),target(ind),dft,cost);
Nlabelled=N

input('next>');

[N,d]=size(data)
ind=[1:N];
fragind=[1:d];

ind=find(max(data,[],2)<=5);
fragind=[2:4];
data=alist(ind,fragind);
target=hlist(ind);
[N,d]=size(data)
disp('before data2lookup');
[lookup,dft,conflict]=data2lookup(data,target);
nconflict=length(conflict)
disp('before catclass');
[lookup,pot,confusion]=catclass(data,target,dft,cost);
Nlabelled=N

input('next>');

  % signature class prediction +++++ against -----
  load bondRT.mat
  ind=(max(alist,[],2)<=20); % remove uncertain entries
  alist=alist(ind,:);
  blist=blist(ind,:);
  N=size(alist,1)
  bclass=ones(N,1);
  signature(1,:)='-----';
  plussign='++++++++++++++++';plussign=plussign';
  ind=find(blist(:,2)>0); % bclass 2 (left 1)
  bclass(ind)=bclass(ind)+1;
  signature(2,:)='+----';
  ind=find(blist(:,3)>0); % bclass 3-4 (left 1,2)
  bclass(ind)=bclass(ind)+2;
  signature([3:4],:)=signature([1:2],:);
  signature([3:4],2)=plussign([1:2]);
  ind=find(blist(:,4)>320); % bclass 5-8 (left 1-4)
  bclass(ind)=bclass(ind)+4;
  signature([5:8],:)=signature([1:4],:);
  signature([5:8],3)=plussign([1:4]);
  ind=find(blist(:,5)>320); % bclass 9-16 (left 1-8)
  bclass(ind)=bclass(ind)+8;
  signature([9:16],:)=signature([1:8],:);
  signature([9:16],4)=plussign([1:8]);
  ind=find(blist(:,6)>320); % bclass 17-32 (left 1-16)
  bclass(ind)=bclass(ind)+16;
  signature([17:32],:)=signature([1:16],:);
  signature([17:32],5)=plussign([1:16]);

  disp('before data2lookup');
  [lookup,dft,conflict]=data2lookup(alist,bclass,0);
  nconflict=length(conflict)

  data1=alist(bclass==1,:);
  data2=alist(bclass==32,:);
  N1=size(data1,1);N2=size(data2,1);
  data=[data1;data2];
  target=[ones(N1,1);1+ones(N2,1)];
  cost=1;
  [N,d]=size(data)

  disp('before data2lookup');
  [lookup1,dft1,conflict]=data2lookup(data,target,0);
  nconflict=length(conflict)

  disp('before catclass');
  [lookup,pot,confusion]=catclass(data,target,dft,cost);
  Nlabelled=N



disp('*** if this works, repeat with sprot distribution ***');








