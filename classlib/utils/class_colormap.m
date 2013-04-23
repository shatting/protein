function [ colormap ] = class_colormap( ncl, method )
%CLASS_COLORMAP(ncl, method=1,2) gets a colormap for ncl classes 

ncl = double(ncl);

if (nargin == 2 && method == 1)
    colormap = jet(ncl);
elseif nargin == 2 && method == 2
    ncl = ncl + 1; % omit black

    s=fix(ncl^(1/3));
    colormap=zeros((s+1)^3,3);
    len=0;
    for a=0:s,
      for b=0:s,
        for c=0:s, 
          len=len+1;
          colormap(len,:)=[a,b,c];
        end;
      end;
    end;
    [int,perm]=sort(sum(colormap(1:len,:),2));
    colormap=colormap(perm(ncl:-1:1),:)/s; % dark part of colors

    colormap = colormap(1:end-1,:); % omit black
else
    colormap = hsv(ncl);
end