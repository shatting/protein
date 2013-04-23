%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% specview.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [sizes,x,y]=specview(data,row,col,rowname,colname)
% spectral line view of a collection of 1D data 
% stratified by two categorical variables
% (many rows are meaningful, but not many columns)
%
% data(l)	l-th data point
% row(l)	row containing l-th data point
% col(l)	column containing l-th data point
%		(or 1 for single column)
% rowname(r,:)  name of row r
% colname(c,:)  name of column c
%
% sizes(r,c)	number of data points in cell (r,c)
% x,y		vectors for reconstructing the plot via
%                  for i=1:size(x,1),
%                    plot(x(i,:),y(i,:),'k-');
%                  end;
%
function [sizes,x,y]=specview(data,row,col,rowname,colname)


smax=20;   % maximal number of spectral lines per spectrum
smax=10;   % maximal number of spectral lines per spectrum




data=data(:)';
n=size(data,2);
if length(row(:))~=n | length(col(:))~=n,
  size_data=size(data)
  size_row=size(row)
  size_col=size(col)
  error('sizes do not match');
end;
x=zeros(2,n);
y=zeros(2,n);

nr=max(row(:));
nc=max(col(:));
if nc==1, col=1+0*row; end;
sizes=zeros(nr,nc);
nlines=0;
for r=1:nr,r
  for c=1:nc,
    ind=find(row==r & col==c);
    n=length(ind);
    sizes(r,c)=n;
    val=sort(data(ind));
    if n==0, 
      % empty data set
      s=0;
    elseif n<smax,
      % extrema only
      s=2;
      val=val([1,n]);
    else
      % smax spectral lines
      s=smax;
      val=val(1+fix((n-0.99999)/s)*[0:s]);
      s=s+1;
    end;
    ind=nlines+[1:s];
    x(1,ind)=c+zeros(1,s);
    x(2,ind)=val;
    y(1,ind)=r+zeros(1,s);
    nlines=nlines+s;
  end;
end;
x=x(:,1:nlines);
y=y(:,1:nlines);

cla;hold on;axis ij    % count rows from top to bottom
% draw boundaries
if nargin>3, 
  xlim=[1,nc+1.5];
  for r=1:nr,text(nc+1.1,r,rowname(r,:));end;
else         
  xlim=[1,nc+1];
end;
if nargin>4, 
  ylim=[-0.2,nr+0.5];
  for c=1:nc,
    h=text(c+0.5,0,colname(c,:));
    set(h,'hor','center');
  end;
else         
  ylim=[0,nr+0.5];
end;
if nc==1,
  set(gca,'ylim',ylim);
else
  % scale values to [0.1,0.9]
  x(2,:)=x(2,:)-min(x(2,:));
  x(2,:)=-0.4+0.8/max(x(2,:))*x(2,:); 
  set(gca,'xlim',xlim,'ylim',ylim);
  for c=1:nc+1,
    plot([c c],ylim,'k-');
  end;
end;
% shift by column index 
x(1,:)=0.5+x(1,:)+x(2,:);
x(2,:)=x(1,:);
% draw spectra
y(2,:)=y(1,:)+0.4;
y(1,:)=y(1,:)-0.4;
plot(x,y,'k');
set(gca,'visible','off')
