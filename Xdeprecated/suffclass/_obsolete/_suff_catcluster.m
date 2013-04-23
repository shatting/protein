%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% catcluster.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lookup,pairfreq,pot]=suff_catcluster(suff,pot)
% function [lookup,pairfreq,pot]=catcluster(dft,pot);
% clustering of categorical data defined by the data frequency table dft
% using an initial potential pot (or the desired number of clusters)
% based on alternating cat2class.m and cat2freq.m
%
% dft(y1,y2,...)        number of occurences of item (y1,y2,...)
% pot(yi,yk,g,ik)       potential for y with y_i=yi,y_k=yk in group g 
%                       ik=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 etc.
%                       i =1 1 2 1 2 3 1 2 3  4  1  2  3  4  5 
%                       k =2 3 3 4 4 4 5 5 5  5  6  6  6  6  6 
% pot=ng                only specifies the desired number of clusters
%                       (maximum 255)
% lookup                group dictionary; see cat2group.m
% pairfreq(yi,yk,g,ik)  frequency of y with y_i=yi,y_k=yk in group g
% pot(yi,yk,g,ik)       revised potential 
%

% suff       sufficient statistic
%            .ncl      number of classes
%            .typ      specifies type of variables
%            .ncat     number of categorical variables
%            .nnum     number of numerical variables
%            .nord     number of ordinal variables
%            .num      list of numerical features
%            .cat      list of categorical features
%            .ord      sublist of ordinal features
%                      ord(cat) is the list of ordinal features
%            .moment   symmetric moment matrix for numerical features
%                      moment(:,:,cl)= sum [x 1]^T*[x 1] 
%                      over all x=data(l,num) with class(l)=cl
%                      moment(end,end,:) is the class size
%                                        even if only categorical data
%            .count    count for categorical features
%                      count(i,a,cl) = number of l with 
%                                      class(l)=cl, data(l,cat(i))=a

ymax=suff.typ(suff.cat);

d=length(ymax);
if d<2, error('data must have more than two components'); end;
ntot=prod(ymax);
yy=max(ymax);

if nargin<3, report=inf; end;

% check potential
potsize=size(pot);
if max(potsize)==1,
  ng=pot;
  if ng>255, 
    ng,
    warning('ng reduced to 255');
    ng=255; 
  end;
  % create random potential
  pot=rand(yy,yy,ng,d*(d-1)/2);
else
  ng=size(pot,3);
  if ng>255, 
    ng,
    error('at most 255 groups are permitted'); 
  end;
  if potsize(4)~=d*(d-1)/2, 
    d,
    error('sizes do not match');
  end;
end;

% initial dictionary
report=1000; % interval length for reporting showtime
[lookup,potsum]=catlookup(ymax,pot,report);
potential=potsum(:)'*dft(:);
disp([showtime,' after intial classification']);

% iteration
it=0;
exchanged=inf;
while exchanged>0,
  it=it+1;
  % cluster sizes
  ng=double(max(lookup.g(:)));
  on=ones(1,ng);
  nglist=zeros(2,ng);
  for g=1:ng,
    ind=find(lookup.g(:)==g);
    nglist(1,g)=length(ind);
    nglist(2,g)=sum(dft(ind));
  end;
  nglist(:,nglist(1,:)>0)

  % compute group frequencies
  pairfreq=cat2freq(suff.count,lookup); 
  save catcluster.mat pot lookup pairfreq
  disp([' it=',num2str(it),': ',showtime,' after catfreq']);

  % get new potential
  pairftot=sum(pairfreq,3);     % sum over all groups
  pairftot=pairftot(:,:,on,:);  % replicate numgroups times
  ind=(pairftot>0);
  pot(ind)=-log(max(pairfreq(ind),realmin)./pairftot(ind));

  % get new dictionary
  lookupold=lookup;
  potentialold=potential;
  [lookup,potsum]=catlookup(ymax,pot);
  potential=potsum(:)'*dft(:);
  gain=potentialold-potential
  exchanged=sum(lookup.g(:)~=lookupold.g(:))
  disp([' it=',num2str(it),': ',showtime,' after catclass']);
end;
