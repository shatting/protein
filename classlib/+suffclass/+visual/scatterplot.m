%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% scatterplot.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function scatterplot(data,class,classname,compname,ndraw, markersize, legend);
% draw a thinned-out scatterplot with up to ndraw points per class
% for multidimensional data, all pair projections are drawn
% 
% data(l,:)	  l-th data point
% class(l,1)	  class of l-th data point 
%                 but class=1 for single class
% classname(c,:)  name of class c  (or '' for default c)
% compname(k,:)   name of component k (or '' for default k)
% ndraw           maximal number of points per class drawn
% markersize      marker size (1=smallest)
%
function scatterplot(data,class,classname,compname,ndraw, markersize, legend, normalize, nclasses);

if (nargin < 7)
    legend = 0;
end
if (nargin < 8)
    normalize = 1;
end

% transpose if necessary
[N,dim]=size(data);
if N==0, error('empty data set'); end;
if dim==1, error('data must be at least 2-dimensional');end;
if length(class)==1,
  ncl=1;
  class=ones(N,1);
else 
    if (nargin >= 9)
       ncl = nclasses; 
    else
       ncl=double(max(class));
    end
  if N~=length(class),
    if dim==length(class),
      data=data';
      [N,dim]=size(data);
    elseif length(class)==1,
      class=ones(N,1);
    else
      N,dim,size(class)
      error('sizes do not fit');
    end;
  end;
  if N~=size(class,1), 
    class=class'; 
  end;
end;

if nargin<5, ndraw=round(500/ncl); end;
if nargin<6, markersize=7; end;

if isempty(classname),
  if ncl>999, error('too many classes for 3i'); end;
  for g=1:ncl,
    classname(g,:)=sprintf('%3i',g);
  end;
end;
if isempty(compname),
  if dim>999, error('too many components for 3i'); end;
  for k=1:dim,
   compname(k,:)=sprintf('%3i',k);
  end;
end;


% prepare class labels and colormap
type ='o+*sdv^<>ph';
ntype=min(ncl,length(type)); % number of different drawing symbols used
while length(type)<ncl, type=[type type]; end;
colormap = class_colormap(ncl);

% select relevant subset of data
for g=ncl:-1:1,  % so that biggest clouds are most visible
  ind=find(class==g);
  nn=length(ind);
  if nn>ndraw,
    % use random selection of ndraw points
    perm=randperm(nn);
    class(ind(perm(ndraw+1:nn)))=0;
  end;
end;
data=double(data(class>0,:));
class=class(class>0);

% plot data
rsub=dim-1;csub=dim-1;if dim==2,csub=2;end;
%clf;
for k=1:dim,
  % normalize data 
  if (normalize)
    y=data(:,k);y=y-min(y);y=y/(max(y)+realmin); 
    data(:,k)=y;  
  else
    y = data(:,k);
  end
  for i=1:k-1,
     if (legend || 1)
        subplot(rsub,csub,i+(dim-1)*(k-2))
     end
    set(gca,'xlim',[0 1],'ylim',[0 1]);hold on;
    % set(gca,'visible','off');
    %set(gca,'xtick',[],'ytick',[]);
    x=data(:,i); % already normalized
    for g=ncl:-1:1,  % so that biggest clouds are most visible
      ind=find(class==g);
      h=plot(x(ind),y(ind),type(g));
      set(h,'markersize',markersize,'color',colormap(g,:));
    end;
%     if i==1,
%       % add names to columns
%       h=text(-0.05,0.5,compname(k,:));
%       set(h,'hor','center');
%     end;
%     if k==dim,
%       % add names to rows
%       h=text(0.5,-0.1,compname(i,:));
%       set(h,'hor','right');
%     end;
  end;
end;

set(gca,'XTick',[]);
set(gca,'YTick',[]);

% identify colors and corresponding classes
if legend,
    
    subplot(rsub,csub,csub);
    set(gca,'xlim',[-3 1],'ylim',[0.5 ncl+0.5]);hold on;
    set(gca,'visible','off');
    image([1:ncl]')
    set(gcf,'colormap',colormap);
    for g=1:ncl,
      h=plot(0.25,g,type(g));
      set(h,'color',colormap(g,:));
      h=text(-0.25,g,classname(g,:));
      set(h,'hor','left');
    end;
end
  
