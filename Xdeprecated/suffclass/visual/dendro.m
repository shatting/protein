% DENDRO Draw a dendrogram.
%
% function g2c=dendro(weight,merge,cost,D,name);
% create dendrogramm (in the unit square) 
% for a merging scheme created by hierarchical clustering 
% (e.g., by clusterh.m)
%
% weight(1,g)         number of vectors (or sum of weights) in group g
% merge(:,i)=[l;k]    merge in step i group l into group k<l 
% cost(1,i)           cost of merging in step i
% D(f,g)              cost of merging groups f and g
%                     (if no preference, D=1);
% name(g,:)           name of group g (default: group number)
%
% g2c(g,nc)           class of group g for a total of nc clusters
%

function g2c=dendro(weight,merge,cost,D,name);

prt=0;                % printlevel (0=nothing,1,2=more,3=next>)

ng=max(merge(:));
if nargin<5,
  name=[];
elseif ~isempty(name) & size(name,1)~=ng, 
  ng,number_of_names=size(name,1)
  error('namelist does not match number of groups');
end;

% create tree with cluster sizes
lchild=zeros(1,2*ng-1);
rchild=zeros(1,2*ng-1);
current=[1:ng];      % current list of nodes
weight(2*ng-1)=0;    % extends the vector to size 2*ng-1
countg(2*ng-1)=0;    % extends the vector to size 2*ng-1
countg(1:ng)=1+countg(1:ng); % unit weight needed for nleft below
for i=ng+1:2*ng-1,
  % merge in step (i-ng) group l into group k<l 
  l=merge(1,i-ng);
  k=merge(2,i-ng);
  lchild(i)=current(k);
  rchild(i)=current(l);
  weight(i)= weight(current(k))+ weight(current(l));
  countg(i)= countg(current(k))+ countg(current(l));
  current(k)=i;
  current(l)=NaN;
  % current % for debug only
end;
if prt>1,children=[lchild;rchild],countg,end;

% assign intervals of labels to the clusters
lbound(2*ng-1)=1;
ubound(2*ng-1)=ng;
for i=2*ng-1:-1:ng+1,
  left=lchild(i);
  right=rchild(i);
  nleft=countg(left);
  lbound(left)=lbound(i);
  ubound(left)=lbound(i)+nleft-1;
  lbound(right)=lbound(i)+nleft;
  ubound(right)=ubound(i);
  if prt>2,
    i,left,right,nleft,bounds=[lbound;ubound]
    input('next>'); 
  end;  
end;
bound=lbound(1:ng);
group(bound)=1:ng; % defines the group label of a bound
if sum(bound~=ubound(1:ng))>0,
  bounds=[lbound;ubound]
  error('programming error: leaves should have lbound=ubound');
end;
if prt>1,bounds=[lbound;ubound],group,end;

% swap branches to decrease distance of neighbors
if length(D)<ng,
  % no preference, no swap 
  ilist=[];
else             
  ilist=ng+1:2*ng-1;
  % disp('swapping is not yet ok, hence currently not done')
  % ilist=[];
end;
for i=ilist,
  left=lchild(i);ll=lbound(left);ul=ubound(left);
  right=rchild(i);lr=lbound(right);ur=ubound(right);
  if D(group(ul),group(lr))>D(group(ur),group(ll)),
    if prt>0, disp(['swap in step ',num2str(i)]); end;
    lgroup=group(ll:ul);
    rgroup=group(lr:ur);
    ul=ll+(ur-lr);
    lr=ul+1;
    group(ll:ul)=rgroup;
    group(lr:ur)=lgroup;
    aux=left;left=right;right=aux;
    lchild(i)=left;
    rchild(i)=right;
    if prt>2, 
      group,children=[lchild;rchild],
      input('next>');
    end;
  end;
end;
if lbound(2*ng-1)==1 & weight(lchild(i))<weight(rchild(i)),
  % swap whole tree to get thick branch to the left (=bottom)
  if prt>0, disp(['swap whole tree']); end;
  aux=lchild;lchild=rchild;rchild=aux;
end;


% reassign intervals of labels to the clusters
lbound(2*ng-1)=1;
ubound(2*ng-1)=ng;
for i=2*ng-1:-1:ng+1,
  left=lchild(i);
  right=rchild(i);
  nleft=countg(left);
  lbound(left)=lbound(i);
  ubound(left)=lbound(i)+nleft-1;
  lbound(right)=lbound(i)+nleft;
  ubound(right)=ubound(i);
  if prt>2, 
    i,left,right,nleft,bounds=[lbound;ubound]
    input('next>'); 
  end;
end;
bound=lbound(1:ng);
group(bound)=1:ng; % defines the group label of a bound
if sum(bound~=ubound(1:ng))>0,
  bounds=[lbound;ubound]
  error('programming error: leaves should have lbound=ubound');
end;
if prt>1,bounds=[lbound;ubound],group,end;


% create relative size and merging costs
size100=round(100*weight/norm(weight,inf));
mincost=min(cost);
if mincost<0, cost=cost-mincost; end; % thus cost>=0;
sumcost=sum(cost);
cost=cost/sumcost;
cost100=1+fix(100*cost/norm(cost,inf));cost100=min(cost100,100);
cost10=1+fix(10*cost/norm(cost,inf));cost10=min(cost10,10);

% assign coordinates y(i) to groups
dy=zeros(1,ng-1);   % y(i+1)-y(i)
for i=1:ng-1,
  ul=ubound(lchild(ng+i));
  dy(ul)=cost10(i);
end;
if min(dy)<=0,
  dy
  error('programming error: dy should be uniquely filled');
end;
y=zeros(1,2*ng-1);
for i=1:ng-1,
  y(group(i+1))=y(group(i))+dy(i);
end;
% normalize to domain [0,1]
y=y/max(y);

% compute coordinates of cluster nodes
x=zeros(1,2*ng-1);
for i=ng+1:2*ng-1,
  left=lchild(i);
  right=rchild(i);
  x(i)=max(x(left),x(right))+1;
  y(i)=(y(left)+y(right))/2;
end;
% reverse and scale x-values
x=1-(x+1)/(max(x)+2);

% plot the dendrogramm
cla;hold on;
set(gca,'xlim',[0 1],'ylim',[0 1],'visible','off');
for k=1:ng,
  g=group(k);
  sizeinfo=[' n=',num2str(size100(g))];
  if isempty(name), text(x(g),y(g),[sizeinfo,';g=',num2str(g)]);
  else              text(x(g),y(g),[sizeinfo,';g=',name(g,:)]);
  end;
end;
for i=ng+1:2*ng-1,
  left=lchild(i);
  right=rchild(i);
  xi=x(i);xl=x(left);xr=x(right);
  yl=y(left);yr=y(right);
  plot([xr,xi,xi,xl],[yr,yr,yl,yl]);
  text(x(i),y(i),...
       [' n=',num2str(size100(i)),'; c=',num2str(cost100(i-ng))]);
  % i,input('next>'); % for debug only
end;
% full cluster, i=2*ng-1
yi=y(i);
plot([xi,xi-1],[yi,yi]);
text(0,1,['n=relative cluster size';'c=relative merging cost']);

% coordinates= [x;y] % for debug only



if nargout==0, return; end;




% create group to class table
g2c=zeros(ng,ng);
nc=ng;
ccurrent=[1:ng]';
g2c(:,nc)=ccurrent;
for i=1:ng-1,
  % merge in step i group l into group k<l 
  l=merge(1,i);
  k=merge(2,i);
  ind=(ccurrent==l);
  ccurrent(ind)=k+0*ccurrent(ind);
  nc=nc-1;
  cnew=ccurrent;
  cc=0;
  for c=1:ng,
    ind=find(cnew==c);
    if ~isempty(ind),
      cc=cc+1;
      cnew(ind)=cc+0*cnew(ind);
    end;
  end;
  g2c(:,nc)=cnew;
end;

