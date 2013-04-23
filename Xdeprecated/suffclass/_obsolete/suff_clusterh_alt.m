%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% suff_clusterh.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [merge,cost,Dsuff0,Dsuffend]=suff_clusterh(type,data,cl,options)
% [merge,cost,Dsuff0,Dsuffend]=suff_clusterh(suff,data,options)
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
    dprintf('creating suff');
    suff = suffstat(max(cl), type);
    suff = suffstat(suff, data, n, cl);
else
    suff = type; % rename
    options = cl;
end

ng = double(suff.ncl);
 
dprintf('\n%i classes, %i datapoints of dimension %i\n',ng,n,d);

if ng<=1,
  % no split
  merge=[];cost=[];sep=[];
  return;
end;

merge=zeros(2,ng-1);
cost=zeros(1,ng-1);

% get regterm, independent of i, j, cl, ncl
reg=1+suff.nnum;  % regularization factor; required reg>1
maxcat = max(suff.typ(suff.cat));
on=ones(1,maxcat,1);
regterm=(1./suff.typ(suff.cat)')*reg(on);
% set positions correctly
shorts = find(suff.typ(suff.cat)< maxcat);
for s=shorts,
    regterm(s,suff.typ(suff.cat(s))+1:maxcat) = 0; % set regterm zero for never ocurring categories
end

% initialize cost matrix
% get datapot
dprintf('getting initial potential values');
[pot,potmin,freq]=classpot(suff);
[cpred, prob, datapot] = classpot(pot, data, size(data,1), freq, options);

% debug
persistent Dsuff;
%persistent mpotijs;
if 1,    
    Dsuff=[];
    %mpotijs=[];
end

% get those out of struct for speed
count = suff.count;
ncat = suff.cat;
cat = suff.cat;

if (isempty(Dsuff)),
    dprintf('calculating initial cost matrix');
    Dsuff = eye(ng);
    Dsuff(eye(ng)==1) = inf;
    
    %mpotijs = cell(double(ng));
    
    % cost matrix

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
%           slow style  
%             % create merged suff        
%             msuff = suffstat(suff,[i,j]); % i+j -> i
% 
%             % get potential for that
%             [mpot,mpotmin,mfreq]=classpot(msuff);
% 
%             [mcpred, mprob, mdatapot] = classpot(mpot, data, size(data,1), options.cat_ent_weight*mfreq);            
            
            % quick&dirty style
            mpotij = mergepot(i,j); % uses count, ncat, cat, regterm, reg, data from parent funtion
            
            idx = setdiff(1:ng,[i,j]);

            mvpred = sum(min([datapot(:,idx) mpotij],[],2));

            Dsuff(i,j) =  mvpred;
        end
        tim = toc;
    end

    Dsuff = Dsuff + Dsuff'
end

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
  
  % add merged classes in count
  count(:,:,k) = count(:,:,k) + count(:,:,l);
  % set mergee counts to 0
  count(:,:,l) = zeros(d,maxcat);  
  % update class k in datapot (using count) (class l will be ignored from now on)
  datapot(:,k) = mergepot(k,l);
  % remove l from index of available classes
  ind(ind==l)=[];
  
  % update cost matrix, k-th row and column
  if (i==1)
      ts = '';
  else
      numleft = (ng-i+1)*(ng-i)/2;
      tim1 = tim/(ng-i); %last iteration
      ts = sprintf(',\t\teta=%s',showtime(tim1*numleft));
  end  
  dprintf('updating cost matrix for merged class %i%s\n',k,ts);
  
  % for all available classes j except k
  for j=setdiff(ind,k),    
    
    % get data potentials for class [k,j]
    mpotkj = mergepot(k,j);

    % remove data potentials for classes k and j
    idx = setdiff(ind,[k,j]);
    
    % insert new data potentials of merged class kj
    % and calculate minimums for each point, and sum
    mvpred = sum(min([datapot(:,idx) mpotkj],[],2));

    % this is the new cost of merging k and j
    Dsuff(k,j) =  mvpred;
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

    function potcl = mergepot(i,j)
	% uses count, ncat, cat, [regterm, reg, data, n {not changing in loops}] from parent funtion
    % doesnt write to any variable from parent
        freq = sum(count,2);
        freq = squeeze(freq(1,:,:));
        
        ncl = sum(freq>0);

        % get potential: add count_i and count_j
        countcl=sum(count(:,:,[i,j]),3)+regterm;        
        countsum=sum(countcl,2);
        mccat=-log(countcl./countsum(:,on));

        % get reg. freq
        freq=freq+reg;
        countall=n+double(ncl)*reg;
        relfreq=squeeze(freq)/countall;

        % get entropy: add i and j rel. frequencies
        potcl=zeros(1,n)-options.cat_ent_weight*log(relfreq(i) + relfreq(j) + realmin);
        
        % calculate new potential values to class i
        for p=1:ncat,
            potcl=potcl+mccat(p,data(:,cat(p)));
        end;
        potcl = potcl';
    end

end