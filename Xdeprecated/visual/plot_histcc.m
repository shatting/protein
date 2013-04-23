function plot_histcc( data, classes, nbins )
%HISTCC class colored histogram. important: classes should not overlap

n = size(data,1);
ncl = max(classes);

classfreq = zeros(max(classes),1);

for i=1:ncl,
    classfreq(i) = sum(classes == i);
end;

% prepare class labels and colormap
%type = 'o+*sdv^<>ph';
%ntype = min(ncl,length(type)); % number of different drawing symbols used
%while length(type)<ncl, type=[type type]; end;
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


hold on;

%hist(data,nbins);
for c=1:max(classes),
    hist(data(classes==c), nbins);
        
    h = findobj(gca,'Type','patch','-not','Tag','colored');
    set(h,'FaceColor',colormap(i,:));
    set(h,'Tag','colored');
end

hold off;