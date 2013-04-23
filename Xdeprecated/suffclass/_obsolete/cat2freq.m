
% 20^5 items in 7 classes take 45s on hektor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cat2freq.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function pairfreq=cat2freq(dft,lookup);
% function pairfreq=cat2freq(dft,lookup,0);
% create a pair frequency table pairfreq for categorical data 
% defined by data frequency table dft and dictionary lookup
%
% function pairfreq=cat2freq(data,target,1);
% create a pair frequency table pairfreq for categorical data 
% defined by labelled (data,target) pairs
%
% dft(y1,y2,...)        number of occurences of item (y1,y2,...)
% lookup                group dictionary; see cat2group.m
% data(l,:)             l-th data vector
% target(l,1)           l-th target (= class number to be predicted)
%
% pairfreq(yi,yk,g,ik)  frequency of y with y_i=yi,y_k=yk in group g
%
function pairfreq=cat2freq(dft,lookup,opt)

if nargin<3,opt=0;end;

if opt==1,
  % (data, target) format
  % data=dft;target=lookup; is used implicitly
  % ik=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 etc.
  % i =1 1 2 1 2 3 1 2 3  4  1  2  3  4  5 
  % k =2 3 3 4 4 4 5 5 5  5  6  6  6  6  6 
  ng=double(max(lookup));
  ymax=double(max(dft,[],1));
  d=size(dft,2);
  yy=double(max(ymax));
  pairfreq=zeros(yy,yy,ng,d*(d-1)/2);
  for g=1:ng,
    disp(['group ',num2str(g),'/',num2str(ng)]);
    indg=(lookup==g);
    for k=2:d,
      for yk=1:ymax(k),
        indk=(indg & dft(:,k)==yk);
        ik=(k-1)*(k-2)/2;
        for i=1:k-1,
          ik=ik+1;
          % disp([ik,i,k])
          for yi=1:ymax(i),
            pairfreq(yi,yk,g,ik)=sum(indk & dft(:,i)==yi); 
          end;
        end;
      end;
    end;
  end;
else
  % (dft,lookup) format
  ymax=size(dft);
  d=length(ymax);
  if d<2, error('data must have more than two components'); end;
  ng=double(max(lookup.g(:)));
  yy=double(max(ymax));
  pairfreq=zeros(yy,yy,ng,d*(d-1)/2);
  nitems=prod(ymax);
  index=reshape(1:nitems,ymax);
  % ik=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 etc.
  % i =1 1 2 1 2 3 1 2 3  4  1  2  3  4  5 
  % k =2 3 3 4 4 4 5 5 5  5  6  6  6  6  6 
  ik=0;
  for k=2:d,
    for i=1:k-1,
      ik=ik+1; % pair position
      % disp([ik,i,k])
      % sort pair to the front
      perm=[i,k,1:i-1,i+1:k-1,k+1:d];
      nik=nitems/(ymax(i)*ymax(k));
      indik=reshape(permute(index,perm),ymax(i),ymax(k),nik);
      % find indices for given pairs
      for yi=1:ymax(i),
        for yk=1:ymax(k),
          indy=indik(yi,yk,:);
          freqy=dft(indy);
          cy=lookup.g(indy);
          for g=1:ng,
            pairfreq(yi,yk,g,ik)=sum(freqy(cy==g));
          end;
        end;
      end;
      % disp([showtime,' after pairs i=',num2str(i),' k=',num2str(k)]);
    end;
  end;
end;

