

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% catlookup_test.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test catlookup.m 



ymax=[5 4 6 7]; % must have length 4
ng=3;
yy=max(ymax);
pot=rand(yy,yy,ng,6);
lookup=catlookup(ymax,pot);
Vlist=zeros(ymax(1),ymax(2),ng);

% pot(yi,yk,g,ik)         potential for y with y_i=yi,y_k=yk in class g 
%                         ik=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 etc.
%                         i =1 1 2 1 2 3 1 2 3  4  1  2  3  4  5 
%                         k =2 3 3 4 4 4 5 5 5  5  6  6  6  6  6 
i=0;
for y4=1:ymax(4),for y3=1:ymax(3),
  y34=[y3,y4]
  for y2=1:ymax(2),for y1=1:ymax(1),
    y=[y1 y2 y3 y4];
    V=pot(y1,y2,:,1)+pot(y1,y3,:,2)+pot(y2,y3,:,3);
    V=V+pot(y1,y4,:,4)+pot(y2,y4,:,5)+pot(y3,y4,:,6);
    V=V(:)';
    Vlist(y1,y2,:)=V;
    [Vmin,g]=min(V);
    if sum(double([lookup.g(y1,y2,y3,y4),cat2group(lookup,y)])-g)>0,
      Vlist
      y,V,g,gl=lookup.g(y1,y2,y3,y4),gc=cat2group(lookup,y)
      disp('wrong group');
      return;
    end;
  end;end;
end;end;



