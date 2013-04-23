
% 20^5 items in 32 classes take 1 min on hektor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% catlookup.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [lookup,potsum]=catlookup(ymax,pot,report);
% create a dictionary lookup for categorical data of shape ymax
% with groups defined by the pair potential pot
%
% the dictionary is used by calling cat2group.m
%
% ymax(k)               number of alternatives in position k
% pot(yi,yk,g,ik)       potential for y with y_i=yi,y_k=yk in group g 
%                       ik=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 etc.
%                       i =1 1 2 1 2 3 1 2 3  4  1  2  3  4  5 
%                       k =2 3 3 4 4 4 5 5 5  5  6  6  6  6  6 
% report                (default inf) interval for showtime
%
% lookup                group dictionary; see cat2group.m
% potsum(y1,y2,...)     potential of item (y1,y2,...) in its group
%
function [lookup,potsum]=catlookup(ymax,pot,report)

d=length(ymax);
if d<2, error('data must have more than two components'); end;
ng=size(pot,3);
if ng>255, ng,error('at most 255 classes are permitted'); end;
if nargin<3, report=inf; end;

% create index table
% ik=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 etc.
% i =1 1 2 1 2 3 1 2 3  4  1  2  3  4  5 
% k =2 3 3 4 4 4 5 5 5  5  6  6  6  6  6 
j=0;
for k=2:d,
  for i=1:k-1,
    j=j+1;
    ik(i,k)=j;ii(j)=i;kk(j)=k;
  end;
end;
ikmax=j;
on1=ones(1,ymax(1));
on2=ones(1,ymax(2));


% accumulate potentials
lookup.g=uint8(zeros(ymax(1),ymax(2),prod(ymax(3:d))));
if nargout==2,
  potsum=zeros(ymax(1),ymax(2),prod(ymax(3:d)));
end;
ind1=1:ymax(1);
ind2=1:ymax(2);
if report<inf, time=toc; end;
y=ones(1,d);
% first two indices are handled loop-free
casmax=prod(ymax(3:d));
for cas=1:casmax,
  V1=zeros(ymax(1),1,ng);
  V2=zeros(ymax(2),1,ng);
  V3=zeros(1,1,ng);
  % pair position 12
  if d > 2,
      V=pot(ind1,ind2,:,1); % V(y1,y2,g)
  else
      V=pot(ind1,ind2,:);
  end
  for k=3:d,
    % pair position 1k
    V1=V1+pot(ind1,y(k),:,ik(1,k)); % independent of position 2
    % pair position 2k,
    V2=V2+pot(ind2,y(k),:,ik(2,k)); % independent of position 1
    for i=3:k-1,
      % pair position ik
      V3=V3+pot(y(i),y(k),:,ik(i,k)); % independent of positions 1,2
    end;
  end;
  V1=V1+V3(on1,:,:);
  V2=permute(V2,[2 1,3]);
  V=V+V1(:,on2,:)+V2(on1,:,:);
  % get group assignment
  [Vmin,lookup.g(:,:,cas)]=min(V,[],3);
  if nargout==3,
    potsum(:,:,ca)=Vmin;
  end;

  % get next y
  for k=3:d, % reshape counts the first index up first
    if y(k)<ymax(k), 
      y(k)=y(k)+1;
      y(3:k-1)=ones(1,k-3);
      break;
    end;
  end;
  if rem(cas,report)==0, 
    disp([showtime(toc-time),' after ',num2str(cas),...
                             ' of ',num2str(casmax)]);
  end;
end;

% create group dictionary
lookup.g=reshape(lookup.g,ymax);

% get index factors 
% such that dft(y1,y2,...)=dft(y*fac-shift)
fac=ones(d,1);
shift=0;
for i=1:d,
  if i<d, 
    fac(i+1)=fac(i)*ymax(i);
    shift=shift+fac(i+1);
  end;
end;
lookup.fac=fac;
lookup.shift=shift;

if nargout==3,
  potsum=reshape(potsum,ymax);
end;
