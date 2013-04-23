%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% specclass.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,y]=specclass(data,class,classname,compname)
% spectral line view of a collection of classified data vectors 
% (many rows are meaningful, but not many columns)
%
% data(l,:)	  l-th data point
% class(l)	  class of l-th data point
% classname(c,:)  name of class c
% compname(k,:)   name of component
%
% x,y		  vectors for reconstructing the plot via
%                    for i=1:size(x,1),
%                      plot(x(i,:),y(i,:),'k-');
%                    end;
%
function [x,y]=specview(data,class,classname,compname)

% normalize data
[N,nc]=size(data);
for k=1:nc,
  datamin=min(data(:,k));
  data(:,k)=data(:,k)-datamin;
  datamax=max(data(:,k));
  data(:,k)=data(:,k)/datamax;
end;

if size(data,1)~=size(class,1), 
  size_data=size(data)
  size_class=size(class)
  error('sizes do not match');
end;


row=class(:,ones(1,nc));
col=[1:nc];col=col(ones(N,1),:);
[sizes,x,y]=specview(data,row,col,classname,compname);

