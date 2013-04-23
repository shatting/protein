%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% catclass.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lookup,pot,confusion]=catclass(data,target,dft,cost,costok);
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
%                       tpred=argmin(cost*confusion(:,g))
%

time0=toc;


% check consistency of input
[N,d]=size(data);
if d<2, error('data must have at least two component'); end;
[NN,t]=size(target);
if t~=1,
  if NN~=1,
    error('currently only single targets are allowed');
  else
    target=target';
    [NN,t]=size(target);
  end;
end;
if NN~=N,
  N,d,NN,t 
  error('sizes do not match');
end;
if max(target)<2, 
  error('at least two target classes are needed');
end;

ymax0=max(data,[],1);
ymin0=min(data,[],1);
if sum(ymin0<=0)>0,
  ymin0
  error('data entries must be positive'); 
end;
if ~(isa(data,'uint8')|isa(data,'uint16')|isa(data,'uint8')...
                      |isa(data,'uint16')),
  if sum(rem(data(:),1)>0)>0,
    error('data entries must be integers'); 
  end;   
end;   

if nargin<3,dft=1;end;
if length(dft)==1,
  % default: uniform distribution
  ymax=double(ymax0);
  dft=ones(ymax);
else
  ymax=size(dft);    % shape of categorical data
  dd=length(ymax);
  if dd~=d,
    N,d,ymax
    error('sizes do not match'); 
  end;
  if sum(ymax0>ymax)>0,
    ymax0,ymax
    error('data mismatch - ymax too small'); 
  end;
end;

nt=double(max(target));   % number of classes to be predicted
ntot=prod(ymax);          % total number of distinct items
yy=max(ymax);             % maximal categorical entry
minsize=fix(sqrt(N))      % minimal size for new groups
Nunl=sum(dft(:));         % number of unlabelled items

if nargin<4,cost=1;end;
if length(cost)==1,
  % default: cost 0 for agreement, 1 for confusion
  cost=1-eye(nt);
end;
report=1000; % interval length for reporting in first call to catlookup 

if nargin<5, costok=0; end;


% initial grouping information
%
% since targets are classes, initial groups = targets (up to labeling); 
% otherwise a clustering is needed and corresponding changes later 
% in step 4(different confusion matrix) 
% and step 6 (group against previous grouping)
%
ng=double(max(target)); % number of classes in target
nglist=zeros(1,ng);
for g=1:ng,
  ind=find(target==g);
  nglist(g)=length(ind);
end;
[sizes,perm]=sort(-nglist);
[aux,ng]=max(find(sizes<0));
ng
t2g(perm(1:ng))=[1:ng];     % target to group table
grouping=t2g(target)'; 
nglist=zeros(1,ng);
for g=1:ng,
  ind=find(grouping==g);
  nglist(g)=length(ind);
end;
on=ones(1,ng);
disp('target class sizes (of labelled entries):');
disp(nglist) 


it=1;itstr=[' it=',num2str(it),': '];
time1=-1;
time=toc;
mcbest=inf;
itbest=0;
stalled=0;

while 1,

  % 1. recluster the available targets, starting with current grouping
  %    this step is here empty since targets are the class numbers


  % 2. get pair frequency table, potentia, and dictionary
  %    from data labelled with current grouping
  gind=(grouping>0); % ignore unlabelled items in data
  pairfreq_lab=cat2freq(data(gind,:),grouping(gind),1);
  disp([itstr,showtime(toc-time),' pairfreq_lab only']);time=toc;

  pot=freq2pot(pairfreq_lab);
  [lookup,potsum]=catlookup(ymax,pot,report);
  report=inf; % don't report in later calls to catlookup
  disp([itstr,showtime(toc-time),' raw dictionary only']);time=toc;


  % 3. get confusion matrix for raw grouping
  pairfreq_unl=cat2freq(dft,lookup); 
  ng=max(grouping);
  if max(lookup.g)<ng,
    % last class is not recognized at all
    pairfreq_unl(1,1,ng,1)=0;
  end;
  disp([itstr,showtime(toc-time),' pairfreq_unl only']);time=toc;

  [gpred,confusion,g2t,mc]=cat2group(lookup,data,target,cost);
  disp([itstr,showtime(toc-time),' raw confusion only']);time=toc;
  mccost_lab=sum(mc);


  % 4. correct potential using unlabelled frequencies
  %    several weights are tried to get maximal generalization
  percent=0;
  for percent=[]%[50,25,10,5,2,1], % relative weight of unlabelled data
    pfw=(percent/100)*N/Nunl;  % weighing factor
    pot=freq2pot(pairfreq_lab+pfw*pairfreq_unl);
    lookup=catlookup(ymax,pot);
    disp([itstr,showtime(toc-time),' new dictionary only']);time=toc;

    % get confusion matrix for new grouping
    [gpred,confusion,g2t,mc]=cat2group(lookup,data,target,cost); 
    disp([itstr,showtime(toc-time),' adapted confusion only']);time=toc;
    mccost=sum(mc);

    if mccost<=max(costok,1.05*mccost_lab+0.01*N), break; end;
  end;
  if isempty(percent), 
    disp('no correction applied -change percent if wanted');
    mccost=mccost_lab; 
    percent=0;
  end;
  disp('rows=target classes, columns=grouping by potential');
  disp(confusion);
  percent_costlab_costnew=[percent mccost_lab mccost]

  if mccost<mcbest,
    mcbest=mccost;
    itbest=it;
    stalled=0;
    % save best classifier
    variables='grouping pot lookup gpred confusion g2t mccost it';
    eval(['save catclass.mat ',variables]);
  elseif mccost==inf,
    error('bad cost');
  else
    stalled=stalled+1;
  end;


  % 5. stopping test
  if time1>=0, 
    disp([itstr,showtime(toc-time1),' per iteration']);
    disp([itstr,showtime(toc-time0),' since calling catclass.m']);
  end;
  if mccost<=costok, break; end;
  if stalled>=5, disp('5 times stalled - break');break; end;
  disp('---------------------------------------------------');
  cont=input('continue classification process? (return=yes,0=no)>');
  if ~isempty(cont), break; end;
  it=it+1;itstr=[' it=',num2str(it),': '];
  time=toc;
  time1=toc;
 

  % 6. refine grouping by intersection with gpred
  ngmax=255; % because of uint8 restriction
  conf=getfreq([grouping(gind),gpred(gind)]);
  disp('rows=current grouping, columns=current predictions');
  disp(confusion);

  % sort intersections by size and remove small ones
  ng=size(conf,1);
  [cmax,cind]=max(conf,[],2);
  for g=1:ng,
    % increase largest entry of row to avoid missing a class
    conf(g,cind(g))=inf;
  end;
  [csort,ind]=sort(-conf(:));
  if length(conf(:))>ngmax, conf(ind(ngmax+1:end))=0; end;
  conf(conf<minsize)=0;

  % create new grouping
  oldgrouping=grouping;
  grouping=zeros(N,1);
  nglist=zeros(1,sum(conf(:)>0));
  gg=0;
  for g=1:ng,
    ind=(oldgrouping==g);
    if cmax(g)<minsize,
      % class must be retained (good splits are too small)
      gg=gg+1;
      grouping(ind)=gg;
      nglist(gg)=sum(ind);
    else
      % class will be grouped, minorities are ignored
      for h=find(conf(g,:)>0),
        gg=gg+1;
        indh=(ind & gpred==h);
        grouping(indh)=gg;
        nglist(gg)=sum(indh);
      end;
    end;
  end;
  disp('group sizes (of labelled entries):');
  disp(nglist) 

  % display change of grouping

  ungrouped0=sum(oldgrouping==0);
  grouping0=oldgrouping;
  if ungrouped0>0,
    disp([num2str(ungrouped0),' ungrouped entries in last row']);
    grouping0(oldgrouping==0)=ng+1;
  end;
  ungrouped1=sum(grouping==0);
  grouping1=grouping;
  if ungrouped1>0,
    disp([num2str(ungrouped1),' ungrouped entries in last column']);
    grouping1(grouping==0)=gg+1;
  end;
  confusion=getfreq([grouping0,grouping1]);
  disp('rows=old grouping, columns=new grouping');
  disp(confusion);
  %disp([itstr,showtime(toc-time),' regrouping only']);time=toc;

  if sum(oldgrouping~=grouping)==0,
    % no change in grouping
    break
  end;
end;


% restore best classifier found
disp(['best classifier found in iteration ',num2str(itbest)]);
disp('saved in catclass.mat - move file to preserve it!');
load catclass.mat
mccost
 
disp(['total time: ',showtime(toc-time0)]);
