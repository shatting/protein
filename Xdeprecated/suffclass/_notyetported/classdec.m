

% priors not yet incorporated
% untested on ordinals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% classdec.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dectree,pot,pmin,cpred]=classdec(typ,data,N,in,out,options);
%          create decision tree data of type typ
%          using data(1:N,in) to predict data(1:N,out) 
%          rounding data to integers is assumed to be be harmless!! 
%          The potential pot (see classpot.m for details) allows one
%          to translate class information into data information
% function cpred=classdec(dectree,data,N);
%          classify incomplete data(1:N,:) using the decision tree
%          [this should not be used without preliminary pruning
%           but preferably postprocess with a pure potential method]
% function classdec(dectree);
%          print information about decision tree
%
% typ(i)     type of feature i
%            c>0: categorical variable with c levels
%            c<0: ordinal variable (meaningless mean) with |c| levels
%            c=0: numerical variable (meaningful mean)
% data(l,:)  l-th feature vector 
% N          number of vectors to be used
% in         list of variables on which decisions are based
% out        list of variables to be predicted
%            (may overlap; in=out for unsupervised classification)
% options    options influencing the tree construction
%            options(1)=opt defines recipe in potential training 
%	       opt = 1: C_i = C for all classes i
%	       opt = 2: C_i = sigma(i)*I (I unity matrix)
%	       opt = 3: C_i = diag(sigma(i,:))
%	       opt = 4: general C_i
%            options(2)=clsiz: minimal class size
%            options(3)=prc: percentage of data declared outliers
%            options(4)=prt: printlevel for classdec
%                       1: list decisions and times 
%                       2: also list loop through classes 
%            options(5)=cbysort: handling of categorical questions
%                       1: sort by min potential
%                       0: use simple clustering
%   
% cpred(l,1) class assignment of l-th feature vector
% dectree    information defining the decision tree
%            ncl=dectree(1)         number of classes 
%            e(node) entry point for node (e=2 at root)
%            i=dectree(e); status variable for node with entry point e
%            if i=0: leaf
%              cpred=dectree(e+1);  class assignment
%            if i>0: split by numerical or ordinal variable i 
%              e1=dectree(e+2);     entry point of child1
%              e2=dectree(e+3);     entry point of child2
%              sep=dectree(e+4);    separator for assignment to children
%                                   data(i)<sep: child1, else child2
%            for i<0: split by categorical, not ordinal variable -i
%              e1=dectree(e+1);     entry point of child1
%              e2=dectree(e+2);     entry point of child2
%              ee=dectree(e+3);     auxiliary index 
%              sep=dectree(e+4:ee); separator set for assignment
%                                   data(i) in sep: child1, else child2
% pot        information defining the potential for classes
%            (see classpot.m for details)
% pmin       potential minima at nodes; for pruning
%        
%
function [dectree,pot,pmin,cpred]=classdec(typ,data,N,in,out,options);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==1,
  % function classdec(dectree);
  % print information about decision tree
  dectree=typ; % rename input

 if nargout>0, error('no output for single input argument'); end;

  ncl=dectree(1);
ncl
  entry=2;
  leaf=0;
  while length(entry)>0,
    e=entry(1);entry(1)=[];   % entry point of current node
    i=dectree(e);
disp(' ');e,i
    if i==0,
      % leaf; assign leaf number (negative)
      leaf=leaf+1;
leaf
    elseif i>0,
      % split by numerical or ordinal variable i 
      e1=dectree(e+1);     % entry point of child1
      e2=dectree(e+2);     % entry point of child2
      entry=[entry e1 e2];
      sep=dectree(e+3);    % separator for assignment to children
                           % data(i)<sep: child1, else child2
entry,sep
    else % i<0
      % split by categorical, not ordinal variable -i
      i=-i;
      e1=dectree(e+1);     % entry point of child1
      e2=dectree(e+2);     % entry point of child2
      entry=[entry e1 e2];
      ee=dectree(e+3);     % auxiliary index 
      sep=dectree(e+4:ee); % separator set for assignment
                           % data(i) in sep: child1, else child2
entry,sep
    end;
  end;
dectree 
  return;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==3,
  % function cpred=classdec(dectree,data,N);
  % classify data(1:N,:) using the decision tree
  dectree=typ; % rename input

  ncl=dectree(1);
  entry=2;
  leaf=0;
  node=2+zeros(N,1);
  prob=zeros(N,ncl);
  while length(entry)>0,
    e=entry(1);entry(1)=[];   % entry point of current node
    ind=find(node==e);        % labels of vectors at current node
    nind=size(ind,1);         % number of vectors at current node
    i=dectree(e);
    if i==0,
      % leaf; assign leaf number (negative)
      leaf=leaf+1;
      node(ind)=zeros(nind,1)-leaf;
    elseif i>0,
      % split by numerical or ordinal variable i 
      e1=dectree(e+1);     % entry point of child1
      e2=dectree(e+2);     % entry point of child2
      entry=[entry e1 e2];
      sep=dectree(e+3);    % separator for assignment to children
                           % data(i)<sep: child1, else child2
      ind1=(data(ind,i)<sep);
      node(ind(ind1))=e1+zeros(sum(ind1),1);
      node(ind(~ind1))=e2+zeros(sum(~ind1),1);
    else % i<0
      % split by categorical, not ordinal variable -i
      i=-i;
      e1=dectree(e+1);     % entry point of child1
      e2=dectree(e+2);     % entry point of child2
      entry=[entry e1 e2];
      ee=dectree(e+3);     % auxiliary index 
      sep=dectree(e+4:ee); % separator set for assignment
                           % data(i) in sep: child1, else child2
      e12=e2+zeros(nind,1);e12(sep)=e1+zeros(ee-e-3,1);
      node(ind)=e12(data(ind,i));
    end;
  end;
  % assign classes (into output variable = dectree)
  dectree=-node;
  return;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dectree,pot,pmin,cpred]=classdec(typ,data,N,in,out,options);
% create decision tree for data of type typ
% using data(1:N,in) to predict data(1:N,out) 
data(N+1:end,:)=[];
opt=options(1);    % defines recipe in potential training
                   % opt = 1: C_i = C for all classes i
 	           % opt = 2: C_i = sigma(i)*I (I unity matrix)
	           % opt = 3: C_i = diag(sigma(i,:))
	           % opt = 4: general C_i
clsiz=options(2);  % minimal class size
prc=options(3);    % percentage of data declared outliers
prt=options(4);    % printlevel
cbysort=options(5);% 1: treat categoricals, sorting by potential minimum
                   % 0: treat categoricals using simple clustering
if prc>0, error(['prc>0 not yet inplemented']); end;
if sum(typ(out)~=0)>0,
  % categorical output variables
  disp('clustering cost vector not yet worked out')
  cbysort=1;
end;


cpred=zeros(N,1);
on=ones(N,1);
nvar=length(typ);
ni=length(in);
nodelist=ones(N,1);
node=0;
lastnode=1;
leaf=0;
maxcl=fix(N/clsiz);                  % maximal number of classes
entry=zeros(maxcl,1);                % list of entry points
sizelist=zeros(1,maxcl);             % list of class sizes
e=2;                                 % current entry point
while node<lastnode,
  node=node+1;
  nodeind=find(nodelist==node);      % labels of vectors at current node
  ncurr=size(nodeind,1);             % number of vectors at current node
  splitlim=max(clsiz,0.2*ncurr);     % current minimal splitsize
  if node==1, 
    sizelist(1)=ncurr;
  else
    if prt,
      tosplit=sizelist(node:lastnode);
      classified=N-sum(tosplit);
      fprintf(1,'\n------------------------------------------------\n');
      fprintf(1,'to split: '); 
      disp(tosplit);
      if classified>0,
        fprintf(1,'classified: ');
        disp(classified);
      end;
      disp('------------------------------------------------');
    end;
    if sizelist(node)~=ncurr,
      node,ncurr,sizelist
      error('mismatch of sizes');
    end;
  end;
  if prt, 
    disp(' ');
    disp(['node ',num2str(node)]);
    disp(['class size ',num2str(ncurr)]); 
  end;
  if node>1 & ncurr<clsiz, error('bad class size'); end;
  % find decision variable
  ipick=0;
  gainpick=inf;
  perm=[];
  val=[];
  if ncurr>=2*clsiz, iilist=[1:ni];  % decisions on input variables only
  else              iilist=[];       % no decision; leaf
  end;
  if prt, 
    isleaf=isempty(iilist);
    if isleaf, disp('leaf (class too small)');
    else       disp([num2str(length(iilist)),' decision variables']); 
    end;
  end;
  for ii=iilist,
    i=in(ii);
    clist=data(nodeind,i);
    % create cumulative sufficient statistic
    if typ(i)>0,
      % categorical decision
      ncl=typ(i);
      suff=suffstat(ncl,typ(out)); % reset sufficient statistic
      suff=suffstat(suff,data(nodeind,out),ncurr,clist);
      if cbysort,
        % sort categorical variables according to potential minimum
        [cpot,cmin]=classpot(suff,opt);
        [potsort,perm]=sort(cmin);
        lmoment=cumsum(suff.moment(:,:,perm),3);
        lcount=cumsum(suff.count(:,:,perm),3); 
      else
        % total sufficient statistic and potential minimum
        moment=sum(suff.moment,3);
        var=diag(moment);
        count=sum(suff.count,3);
        suffall=suff;
        suffall.ncl=1;
        suffall.moment=moment;
        suffall.count=count;
        [tpot,tmin,tfreq]=classpot(suffall,opt);
        pmin(node)=tmin-log(tfreq);

        % get Cholesky factor for scale
        mucurr=moment(1:end-1,end)/ncurr;
        cov=moment(1:end-1,1:end-1)/ncurr-mucurr*mucurr';
        cov=cov+diag(var(1:end-1)/ncurr^2); % regularization
        % take care of null rows
        normcov=norm(cov(:),inf);
        nullrows=find(diag(cov)<=0)';
        for row=nullrows, 
          cov(row,row)=1e-8*normcov; 
        end;
        [R,p]=chol(cov);
        if p>0,
          nullrows,cov,eig(cov)
          error('insufficient regularization before cluster');
          d=size(cov,1);R(d,d)=1;
          R(p:d,p:d)=eye(d+1-p);
        end;

        % try split by clustering
        countx=suff.moment(end,end,:);
        countx=countx(:)';
        % numerical target variables only
        sumx=permute(suff.moment(1:end-1,end,:),[1 3 2]);
        [merge,cost,sep]=clusterh(countx,R'\sumx,splitlim);
        lmoment=sum(suff.moment(:,:,sep),3);
        lcount=sum(suff.count(:,:,sep),3); 
      end;
    else % typ(i)<=0, 
      % numerical or ordinal decision
      val=rsort(clist);
      ncl=length(val);
      suff=suffstat(1,typ(out)); % reset sufficient statistic
      for cl=1:ncl,
        if prt>1, disp([showtime,'   cl=',num2str(cl)]); end;
        nodecl=nodeind(clist==val(cl));
        nncl=length(nodecl);
        suff=suffstat(suff,data(nncl,out),nncl,ones(nncl,1));
        lmoment(:,:,cl)=suff.moment(:,:,1);
        lcount(:,:,cl)=suff.count(:,:,1);
      end;
      lmoment=lmoment(:,:,1:ncl); 
      if isempty(lcount), lcount=zeros(0,0,ncl);
      else                lcount=lcount(:,:,1:ncl); 
      end;
    end; % of create cumulative sufficient statistic

    if typ(i)<=0 | cbysort,
      % total sufficient statistic and potential minimum
      moment=lmoment(:,:,ncl);
      count=lcount(:,:,ncl);
      suff.ncl=1;
      suff.moment=moment;
      suff.count=count;
      [tpot,tmin,tfreq]=classpot(suff,opt);
      pmin(node)=tmin-log(tfreq);

      % sufficient statistic for splits and potential gains
      Vgain=inf+zeros(1,ncl);
      lsiz=lmoment(end,end,:);
      cllist=find(lsiz>=splitlim & lsiz<=ncurr-splitlim);
      nsplit=length(cllist);
      if prt, disp(['try ',num2str(nsplit),' splits']); end;
    else
      % clustering already performed, at most one class tried
      if isempty(sep), 
        cllist=[];
      else             
        cllist=1;perm=sep; % store for permpick
      end;
    end;
    for cl=cllist,
      mom=lmoment(:,:,cl);
      ct=lcount(:,:,cl);
      suff.ncl=1;
      suff.moment=mom;
      suff.count=ct;
      [lpot,lmin,lfreq]=classpot(suff,opt);
      suff.moment=moment-mom;
      suff.count=count-ct;
      [rpot,rmin,rfreq]=classpot(suff,opt);
      Vgain=lmin-log(lfreq)+rmin-log(rfreq)-pmin(node);
      if Vgain<gainpick,
        % save improved choice
        gainpick=Vgain;
        ipick=i;
        clpick=cl;
        permpick=perm;
        valpick=val;
        nlpick=mom(end,end);
        nrpick=ncurr-nlpick;
      end;
    end;
    if prt, disp([showtime,'   ii=',num2str(ii)]); end;
  end; % of for ii

  % split or assign class
  if prt, 
    if ipick>0,
      if typ(ipick)==0, disp('numerical split');
      elseif typ(ipick)>0, disp('categorical split');
      else                 disp('ordinal split');
      end;
      disp(['ipick=',num2str(ipick),' nlpick=',num2str(nlpick),...
            ' nrpick=',num2str(nrpick)]);
    elseif ~isleaf,
      disp('leaf (split failed)');
    end;
  end;
  entry(node)=e;
  if ipick==0,
    % leaf; assign class
    leaf=leaf+1;
    cpred(nodeind)=leaf+cpred(nodeind);
    if leaf==1,
      % create potential
      pot=tpot;
      % reserve space
      pot.ccat(:,:,maxcl)=tpot.ccat;
      pot.mu(:,:,maxcl)=tpot.mu;
      pot.R(:,:,maxcl)=tpot.R;
      pot.shift(:,:,maxcl)=tpot.shift;  
    else
      % extend tree potential by new class 
      pot.ncl=leaf;
      pot.ccat(:,:,leaf)=tpot.ccat;
      pot.mu(:,:,leaf)=tpot.mu;
      pot.R(:,:,leaf)=tpot.R;
      pot.shift(:,:,leaf)=tpot.shift;  
    end;
    dectree(e)=ipick;
    dectree(e+1)=leaf;     % class assignments
    e=e+2;                 % new entry point
  elseif typ(ipick)<=0,
    % split by numerical or ordinal variable 
    ch1=lastnode+1;        % node of child1
    ch2=lastnode+2;        % node of child2
    lastnode=lastnode+2;
    % create separator value
    sepl=valpick(clpick);
    sepr=valpick(clpick+1);
    sep=sepl+(sepr-sepl)/2;
    sepround=1+fix(sep);
    if sepround>sepl & sepround<=sepr, sep=sepround; end;
    sizelist(ch1)=nlpick;
    sizelist(ch2)=nrpick;
    if prt, disp(['sep=',num2str(sep)]); end;
    % assign data to children nodes
    ind1=(data(nodeind,ipick)<sep);
    nodelist(nodeind(ind1))=ch1+zeros(sum(ind1),1);
    nodelist(nodeind(~ind1))=ch2+zeros(sum(~ind1),1);
    dectree(e)=ipick; 
    dectree(e+1)=ch1;    
    dectree(e+2)=ch2;   
    dectree(e+3)=sep;     % separator for assignment to children
                          % data(i)<sep: child1, else child2
    e=e+4;                % new entry point
  else 
    % split by categorical, not ordinal variable 
    ch1=lastnode+1;       % node of child1
    ch2=lastnode+2;       % node of child2
    lastnode=lastnode+2;
    % create separator set
    if cbysort,
      perminv(permpick)=[1:ncl]';
      di=rsort(data(nodeind,ipick)); % actual values of variable
      invdi=perminv(di);
      nlsep=sum(invdi<=clpick);
      nrsep=sum(invdi>clpick);
      if nlsep<=nrsep, 
        sep=permpick(invdi(invdi<=clpick));
        sizelist(ch1)=nlpick;
        sizelist(ch2)=nrpick;
      else             
        sep=permpick(clpick+1:end);
        sizelist(ch1)=nrpick;
        sizelist(ch2)=nlpick;
      end;
    else
      % sep already assigned by clustering
      sep=permpick;
      sizelist(ch1)=nlpick;
      sizelist(ch2)=nrpick;
    end;
    if prt, fprintf(1,'separator ');disp(sep); end;
    % assign data to children nodes
    nsep=length(sep);
    node12=ch2+zeros(ncl,1);
    node12(sep)=ch1+zeros(nsep,1);
    nodelist(nodeind)=node12(data(nodeind,ipick));
    if 0, % for debug only
      nch1=sum(nodelist==ch1)
      nch2=sum(nodelist==ch2)
      if nch1~=nlpick & nch1~=nrpick, error('bad count'); end;
    end;
    dectree(e)=-ipick;
    dectree(e+1)=ch1;     % entry point of child1
    dectree(e+2)=ch2;     % entry point of child2
    ee=e+3+nsep;
    dectree(e+3)=ee;     % auxiliary index 
    dectree(e+4:ee)=sep; % separator set for assignment
                         % data(i) in sep: child1, else child2
    e=ee+1;              % new entry point
  end; % of split or assign class

end; % of while unprocessed node


% assign number of classes
dectree(1)=leaf;

% replace child node numbers by entry points
for node=1:lastnode,
  e=entry(node);
  i=dectree(e);
  if i~=0,
    dectree(e+1:e+2)=entry(dectree(e+1:e+2));
  end;
end;


% remove unused space
if pot.ncat>0,
  pot.ccat=pot.ccat(:,:,1:leaf);
end;
if pot.nnum>0,
  pot.mu=pot.mu(:,:,1:leaf);
  pot.R=pot.R(:,:,1:leaf);
  pot.shift=pot.shift(:,:,1:leaf);  
end;
if prt, disp(' ');disp(' ');
else    lastnode
end;
