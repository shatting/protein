% SUFFSTAT Create, add and remove data to or from a sufficient statistic.
%
%   A sufficient statistic summarizes information about data vectors.
%   Vectors may consist of any number of numerical (real numbers), 
%   categorical (integers with meaningless mean value) or ordinal (integers 
%   with meaningful mean value) features. A suffstat structure is 
%   initialized with a call to
%
%       suff=suffstat(ncl,typ)
%   
%   where ncl, the number of classes desired, is an integer, and typ is a
%   1x(number of features) vector describing the type of each feature:
%   
%       typ(i)>0: categorical variable (meaningless mean) with typ(i) levels
%       typ(i)<0: ordinal variable (meaningful mean) with |typ(i)| levels
%       typ(i)=0: numerical variable
%   
%   Once initialized, the structure has to be passed on to all calls that 
%   should operate on the statistic it represents. For instance, the
%   following call adds information from the vectors contained in data to
%   the statistic suff:
%
%       suff=suffstat(suff,data,N,clist)
%
%   clist is a vector containing class labels for data. Please note, 
%   that any vectors with lables not falling into the range 1:suff.ncl will 
%   be ignored in the call and will not be added to the statistic.
%
% function suff=suffstat(ncl,typ)
%       initialize sufficient statistic for data of arbitrary type 
%       (>=1 feature) and ncl classes
%
% function suff=suffstat(suff,data,N,clist,[],wt);
%       increment sufficient statistic using data(1:N,:) with assigned 
%       classes clist(1:N) and weights wt
%
% function suff=suffstat(suff,data,N,clist,cold,wt,wtold);
%       modify sufficient statistic using data(1:N,:) when class 
%       assignments change from cold(1:N) to clist(1:N) with new weights 
%       wt and old weights wtold
%
% function suff=suffstat(suff,ncl);
%       increase the number of classes to ncl
%
% function suff=suffstat(suff,mergelist);
%       merge classes from mergelist 
%
% INPUTS:
%   ncl           number of classes
%   typ(i)        type of feature i
%                   typ(i)>0: categorical variable (meaningful mean) with 
%                             typ(i) levels
%                   typ(i)<0: ordinal variable (meaningful mean) with 
%                             |typ(i)| levels
%                   typ(i)=0: numerical variable
%   suff          a sufficient statistic structure, generated from any prior
%                 call to suffsat (see outputs section for details on fields)
%   data(l,:)     l-th feature vector in current batch
%   N             number of vectors in current batch
%   clist(l)      class assigned to l-th feature vector
%   wt(l)         weight assigned to l-th feature vector
%   mergelist     
%
% OUTPUTS:
% suff       sufficient statistic
%            .ncl      number of classes
%            .typ      specifies type of features
%            .ncat     number of categorical features
%            .nnum     number of numerical features
%            .nord     number of ordinal features
%            .num      list of numerical feature indices
%            .cat      list of categorical feature indices
%            .ord      sublist of ordinal feature indices
%                      .ord(.cat) is the list of ordinal feature indices
%            .moment   symmetric moment matrix for numerical features
%                      moment(:,:,cl)=sum [x 1]^T*[x 1] over all 
%                      x=data(l,.num) with class(l)=cl. moment(end,end,cl) 
%                      is size of class cl, even if only categorical data
%            .count    count for categorical features
%                      count(i,a,cl) = number of l with 
%                                      class(l)=cl and data(l,cat(i))=a
%
%

function s=suffstat(s,data,N,clist,cold,wt,wtold)
prt=0;  % printlevel

if isstruct(s) && nargin==2,
  if length(data) == 1,
      % function suff=suffstat(suff,ncl);
      % increase the number of classes to ncl
      ncl=data; % rename input data
      if ncl<=s.ncl, 
        error('there are already enough classes');
      end;
      s.ncl=ncl;
      s.moment(1,1,ncl)=0;
      s.count(1,1,ncl)=0;
   else
     % function suff=suffstat(suff,mergelist);
     mergelist = data; % rename
     newcl = min(mergelist);
     s.moment(:,:,newcl) = sum(s.moment(:,:,mergelist),3);     
     % want the classes to stay in place, eg. for clustering
     s.moment(:,:,mergelist(2:end)) = s.moment(:,:,mergelist(2:end))*0;
     s.count(:,:,newcl) = sum(s.count(:,:,mergelist),3);
     s.count(:,:,mergelist(2:end)) = s.count(:,:,mergelist(2:end))*0;       
  end      
  return;
end;

if nargin==2 
  % function suff=suffstat(ncl,type);
  % initialize sufficient statistic 
  % for data of arbitrary type (>1 feature) and ncl classes
  ncl=s;type=double(data); % rename input data

  num=find(type==0);nnum=length(num);
  cat=find(type~=0);ncat=length(cat);
  typc=type(cat); % restrict type to categorical data
  ord=find(typc<0);nord=length(ord);
  moment=zeros(nnum+1,nnum+1,ncl);
  count=zeros(ncat,max(typc),ncl);
  s=struct('ncl',ncl,'typ',type,'nnum',nnum,'ncat',ncat,'nord',nord,...
           'num',num,'cat',cat,'ord',ord,...
           'moment',moment,'count',count);
  return;
end;

ncl = s.ncl;
ncat = s.ncat;
cat = s.cat;
typ = s.typ;
count = s.count;

% function suff=suffstat(suff,data,N,clist,[],wt);
% increment sufficient statistic using data(1:N,:) 
% with assigned classes clist(1:N)
if (nargin >= 6), withweights=1; else withweights=0;end

for cl=1:ncl,
    if prt, disp(['increment class ',num2str(cl)]); end;
    clind=find(clist(1:N)==cl);
    nclind=length(clind);
    for i=1:ncat,
        if prt, disp(['i=',num2str(i)]); end;
        ii=cat(i);
        for a=1:typ(ii),
            aaidx = data(clind,ii)==a;
            if withweights,
                count(i,a,cl)=count(i,a,cl)+sum(wt(clind(aaidx)));
            else
                count(i,a,cl)=count(i,a,cl)+sum(aaidx);
            end
        end;
    end;
    if s.nnum>0,
        if withweights,
            X=[data(clind,s.num), ones(nclind,1)];  % [x_i 1]
            Y=X.*repmat(wt(clind),1,s.nnum+1);      % [x_i*w_i w_i]
            X=X'*Y;
        else
            X=[data(clind,s.num), ones(nclind,1)];
            X=X'*X;
        end
        s.moment(:,:,cl)=s.moment(:,:,cl)+X;        
    else
        if withweights,
            s.moment(:,:,cl)=s.moment(:,:,cl)+sum(wt(clind));
        else
            s.moment(:,:,cl)=s.moment(:,:,cl)+nclind;
        end
    end;
end;
if nargin<=6, s.count = count; return; end;

% function suff=suffstat(suff,data,N,clist,cold,wt,wtold);
% decrement contribution from old class assignment

for cl=1:ncl,
    if prt, disp(['decrement class ',num2str(cl)]); end;
    clind=(cold(1:N)==cl);
    nclind=length(clind);
    for i=1:ncat,
        ii=cat(i);
        for a=1:typ(ii),
            aaidx = data(clind,ii)==a;
            if withweights,
                count(i,a,cl)=count(i,a,cl)-sum(wtold(clind(aaidx)));
            else
                count(i,a,cl)=count(i,a,cl)-sum(aaidx);
            end
        end;
    end;
    if s.nnum>0,
        if withweights,
            X=[data(clind,s.num), ones(nclind,1)];  % [x_i 1]
            Y=X.*repmat(wtold(clind),1,s.nnum+1);      % [x_i*w_i w_i]
            X=X'*Y;
        else
            X=[data(clind,s.num), ones(nclind,1)];
            X=X'*X;
        end
        s.moment(:,:,cl)=s.moment(:,:,cl)-X;
    else
        if withweights,
            s.moment(:,:,cl)=s.moment(:,:,cl)-sum(wtold(clind));
        else
            s.moment(:,:,cl)=s.moment(:,:,cl)-nclind;
        end
    end;
end;

s.count = count;
