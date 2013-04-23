

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% data2lookup.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [lookup,dft,conflict]=data2lookup(data,target,opt);
% create a dictionary lookup and a data frequency table
% for categorical data given in (data,target) form
%
% data(l,:)             l-th data vector
% target(l,1)           l-th target (single class only)
% opt                   1: count conflicting data in dft (default)
%                       0: ignore conflicting data in dft
%                          (count only agreements with first)
%
% dft(y1,y2,...)        number of occurences of item (y1,y2,...)
% lookup                group dictionary; see cat2group.m
% conflict              conflict list: for l=conflict.
%                       data(l,:) has target distinct
%                       from corresponding lookup entry
%
function [lookup,dft,conflict]=data2lookup(data,target,opt);

[N,d]=size(data);
if d<2, error('data must have at least 2 component'); end;
[NN,t]=size(target);
if t~=1,
  if NN~=1,
    error('only single targets are allowed');
  else
    target=target';
    [NN,t]=size(target);
  end;
end;
if NN~=N,
  N,d,NN,t 
  error('sizes do not match');
end;

if nargin<3, opt=1; end;

% get index factors 
% such that dft(y1,y2,...)=dft(y*fac-shift)
fac=ones(d,1);
shift=0;
for i=1:d,
  ymax(i)=double(max(data(:,i)));
  if i<d, 
    fac(i+1)=fac(i)*ymax(i);
    shift=shift+fac(i+1);
  end;
end;
lookup.fac=fac;
lookup.shift=shift;
lookup.g=zeros(ymax);
dft=zeros(ymax);
nconflict=0;
nconfmax=10000;
conflict=zeros(1,nconfmax);

tic
for l=1:N,
  y=double(data(l,:));
  ind=y*fac-shift;
  freq=dft(ind);
  if freq==0,
    lookup.g(ind)=target(l);
    dft(ind)=1;
  else
    if lookup.g(ind)~=target(l),
      % add to conflict list
      nconflict=nconflict+1;
      if nconflict>nconfmax,
        nconfmax=nconfmax+10000;
        conflict(1,nconfmax)=0;
      end;
      conflict(1,nconflict)=l;
      if opt, dft(ind)=freq+1; end;
    else
      dft(ind)=freq+1;
    end;
  end;
  if rem(l,50000)==0, 
    t=toc;tt=toc*(N-l)/l;
    disp([showtime(t),' for ',num2str(l),'; missing ',showtime(tt)]); 
  end;
end;
conflict=conflict(1:nconflict);
