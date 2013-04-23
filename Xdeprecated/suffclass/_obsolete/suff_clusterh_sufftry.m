%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% suff_clusterh.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [merge,cost,Dsuff0,Dsuffend]=suff_clusterh(type,data,cl,options)
% [merge,cost,Dsuff0,Dsuffend]=suff_clusterh(suff,options)
% hierarchical clustering 
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
function [merge,cost,Dsuff0,Dsuffend]=suff_clusterh(type,data,cl,options)


dprintf('----------------- start of suff_clusterh.m ------------------\n');
[n,d] = size(data);

if ~isstruct(type),
    % suff_clusterh(type,data,cl,options)
    dprintf('creating suff');
    suff = suffstat(max(cl), type);
    suff = suffstat(suff, data, n, cl);
else
    % suff_clusterh(suff,options)
    suff = type; % rename
    options = data; % rename
    
    n = sum(suff.moment(end,end,:)); %? true for mixed suffs?
    d = length(suff.typ);
    dcat = length(suff.cat);
end


potind = zeros(size(suff.count(:,:,1)));
if (any(suff.typ(suff.cat) < max(suff.typ(suff.cat)))),
    %calculate count pattern: we dont want the infs that are outside of the
    %category range to mess up the sum(pot.*count)    
    for i=1:dcat,
       potind(i,1:suff.typ(suff.cat(i))) = 1; 
    end   
else
    potind = potind + 1;
end
potind = logical(potind);


ng = double(suff.ncl);
 
dprintf('\n%i classes, %i datapoints of dimension %i\n',ng,n,d);

if ng<=1,
  % no split
  merge=[];cost=[];sep=[];
  return;
end;

merge=zeros(2,ng-1);
cost=zeros(1,ng-1);

% initialize cost matrix
[pot,potmin,freq]=classpot(suff);

% remember class pot values: we dont want to recalc them every time
for i=1:ng,
    classpotvals(i)= datapot(pot.ccat(:,:,i), suff.count(:,:,i), freq(i));
end

% cost matrix
dprintf('initializing cost matrix');
Dsuff = eye(ng);
Dsuff(eye(ng)==1) = inf;
for i=1:ng-1,
    tic;
    if (i==1)
        ts='';
    else
        numleft = (ng-i+1)*(ng-i)/2;
        tim1 = tim/(ng-i);
        ts = sprintf(',\t\teta=%s.',showtime(tim1*numleft));
    end
    dprintf('pairs [%i,%i:%i]%s',i,i+1,ng,ts);
    for j=i+1:ng,            
        msuff = suffstat(suff,[i,j]);
        [mpot,mpotmin,mfreq]=classpot(msuff);
        % dV_c1c2 = V_c1c2 - V_c1 - V_c2
        Dsuff(i,j) = datapot(mpot.ccat(:,:,i), msuff.count(:,:,i), mfreq(i)) - classpotvals(i) - classpotvals(j);
    end
    tim = toc;
end

Dsuff = Dsuff + Dsuff'

% save initial cost matrix
Dsuff0 = Dsuff;

dprintf('starting merges');
ind=1:ng;
for i=1:ng-1,
  tic;
  % find merge information
  [clist,klist]=min(Dsuff(ind,ind),[],1);
  [c,k]=min(clist);
  l=min(find(clist==c));
  k=ind(k);
  l=ind(klist(l));
  if k<l, merge(1,i)=l;merge(2,i)=k;
  else    merge(1,i)=l;merge(2,i)=k; % interchange would cause trouble?
  end;
  cost(i)=c;
  dprintf('%i/%i-th merge: [%i,%i]->%i at a cost of %f',i,ng,k,l,k,c);
  
  % add classes in suff
  suff=suffstat(suff,[k,l]);
  % update class k in classpotvals (class l will be ignored from
  % now on)
  [mpot,mpotmin,mfreq]=classpot(msuff);
  
  classpotvals(k) = datapot(mpot.ccat(:,:,k), msuff.count(:,:,k), mfreq(k));
  % remove l from index of available classes
  ind(ind==l)=[];
  

  if (i==1)
      ts = '';
  else
      numleft = (ng-i+1)*(ng-i)/2;
      tim1 = tim/(ng-i); %last iteration
      ts = sprintf(',\t\teta=%s',showtime(tim1*numleft));
  end  
  dprintf('updating cost matrix for merged class %i%s\n',k,ts);

  % update cost matrix, k-th row and column  
  % for all available classes j except k
  for j=setdiff(ind,k),       
    msuff = suffstat(suff,[k,j]);
    [mpot,mpotmin,mfreq]=classpot(msuff);
	Dsuff(k,j) = datapot(mpot.ccat(:,:,k), msuff.count(:,:,k), mfreq(k)) - classpotvals(k) - classpotvals(j);
  end    

  % copy row to column
  Dsuff(k,:) = Dsuff(:,k)';
  
  tim = toc;
end;

Dsuffend = Dsuff;

% find separating index set
% if nargin<3, splitlim=0; end;
% countall=countx(ind); % total sum of weights
% splitlim2=countall-splitlim; 
% kok=0;
% lok=0;
% splitok=0;
% ind=zeros(1,ng0);
% for i=ng-1:-1:1,
%   l=merge(1,i);
%   k=merge(2,i);
%   countx(k)=countx(k)-countx(l);
%   if ~splitok,
%     kok=(countx(k)>=splitlim & countx(k)<=splitlim2);
%     lok=(countx(l)>=splitlim & countx(l)<=splitlim2);
%     splitok=(kok|lok);
%     if kok, ind(k)=1;kok=countx(k); end;
%     if lok, ind(l)=2;lok=countx(l); end;
%   else,
%     ind(l)=ind(k);
%   end;
% end;
% if max(kok,lok)==0, sep=[];           wt=0;
% elseif kok>=lok,    sep=find(ind==1); wt=kok;
% else                sep=find(ind==2); wt=lok;
% end;

% split=[wt,countall-wt]

    function pot = datapot(cpot,ccount,pr)
        dfreq = sum(ccount(1,:));
        pot = cpot(potind)'*ccount(potind) - log(pr + realmin)*dfreq;
    end

end