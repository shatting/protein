% SUFF_CLUSTERH Hierarchical cluster analysis of (only categorical yet) data.
%
% [merge,cost,Dsuff0,Dsuffend]=suff_clusterh(type, data, cl, options, verbose)
% hierarchical clustering of row vectors contained in data, having features described by type
% and preassigned classes cl.
%
% [merge,cost,Dsuff0,Dsuffend]=suff_clusterh(suff, options, verbose)
% hierarchical clustering of the classes described by suff
% 
% INPUT:
%
% type(i)             type of feature i
% data(l,:)           l-th data vector
% cl(l)               class of l-th data vector
% 
% suff                suffictient statistic describing the data and classes, see suffstat.m 
%
% options             not used as of now
%
% OUTPUT: 
%
% merge(:,i)=[l;k]    merge in step i class l into class k<l 
% cost(1,i)           cost of merging in step i
% Dsuff0(f,g)         initial merge cost matrix
% Dsuffend(f,g)       final cost matrix
%
% See also: suff_classify.

function [merge,cost,Dsuff0,Dsuffend]=clusterh(type,data,cl,options,verbose)

dprintf('----------------- start of suff_clusterh.m ------------------\n');
[n,d] = size(data);

% build suff from data or rename inputs
if ~isstruct(type),
    dprintf('creating suff');
    suff = suffstat(max(cl), type);
    suff = suffstat(suff, data, n, cl);
    if (nargin<5), verbose = 0; end;
else
    suff = type; % rename
    options = data;
    if (nargin<3), verbose = 0; else verbose = cl;end
    n = sum(suff.moment(end,end,:));
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

% get regterm, independent of i, j, cl, ncl
reg=1+suff.nnum;  % regularization factor; required reg>1
maxcat = max(suff.typ(suff.cat));
on=ones(1,maxcat,1);
regterm=(1./suff.typ(suff.cat)')*reg(on);
% set positions correctly
shorts = find(suff.typ(suff.cat)< maxcat);
for s=shorts,
    regterm(s,suff.typ(suff.cat(s))+1:maxcat) = 0; % set regterm zero for never occuring categories
end

% initialize cost matrix
dprintf('getting initial potential values');
[pot,potmin,freq]=classpot(suff);

% get those out of struct for speed
count = suff.count;
dprintf('calculating initial cost matrix');
Dsuff = eye(ng);
Dsuff(eye(ng)==1) = inf;

%mpotijs = cell(double(ng));

% cost matrix

for i=1:ng,
    pots(i) = datapot(pot.ccat(:,:,i),count(:,:,i),freq(i));
end

for i=1:ng-1,
    tic;
    if (verbose)
        if (i==1)
            ts='';
        else
            numleft = (ng-i+1)*(ng-i)/2;
            tim1 = tim/(ng-i);
            ts = sprintf(',\t\teta=%s.',showtime(tim1*numleft));
        end
        dprintf('pairs [%i,%i:%i]%s',i,i+1,ng,ts);
    end
    for j=i+1:ng,            

        mcount = count(:,:,i) + count(:,:,j);
        [mpot, mfreq] = mergepot(i,j);            
        % dV_c1c2 = V_c1c2 - V_c1 - V_c2
        Dsuff(i,j) = datapot(mpot,mcount,mfreq(i)) - pots(i) - pots(j);

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
  if (verbose)
      dprintf('%i/%i-th merge: [%i,%i]->%i at a cost of %f',i,ng,k,l,k,c);
  end
  
  % add merged classes in count
  count(:,:,k) = count(:,:,k) + count(:,:,l);  
  % set mergee counts to 0
  count(:,:,l) = zeros(d,maxcat);  
  % update class k in datapot (using count) (class l will be ignored from now on)
  [mpot, mfreq] = mergepot(k,l);
  pots(k) = datapot(mpot,count(:,:,k),mfreq(k));
  
  % remove l from index of available classes
  ind(ind==l)=[];
  
  % update cost matrix, k-th row and column
  if (verbose)
      if (i==1)
          ts = '';
      else
          numleft = (ng-i+1)*(ng-i)/2;
          tim1 = tim/(ng-i); %last iteration
          ts = sprintf(',\t\teta=%s',showtime(tim1*numleft));
      end  
      dprintf('updating cost matrix for merged class %i%s\n',k,ts);
  end
  
  % for all available classes j except k
  for j=setdiff(ind,k),       
    mcount = count(:,:,j) + count(:,:,k);
    [mpot, mfreq] = mergepot(j,k);            
   Dsuff(j,k) = datapot(mpot,mcount,mfreq(j)) - pots(j) - pots(k);
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

    function [mccat, relfreq] = mergepot(i,j)
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

    end

end