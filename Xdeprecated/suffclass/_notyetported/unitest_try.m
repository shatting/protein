

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% unitest_try.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [teststat,testcor]=unitest(data);
% computes various test statistics for testing uniformity of data 
% 
function [teststat,testcor]=unitest(data);

if 0,
  % find average support underestimation
  % using one unit on each side is ok in the mean
  % thus we compare later the normalized version against i/(N+1) 
  tic
  clear
  ncas=1000;
  nlist=[10000 5000 2000 1000 500 200 100 50 20 10 5 4 3 2];
  fraclist=[1/ncas,0.005,0.01,0.05,0.1,0.5,0.9,0.95,0.99,0.995];
  % fraclist=[1/ncas,0.005,0.01,0.05,0.1,0.5];
  nn=length(nlist);
  frac=round(fraclist*ncas);
  wid=zeros(nn,ncas);
  widk=zeros(1,ncas);
  for k=1:nn,
    n=nlist(k);
    for cas=1:ncas,
      data=rand(n,1);
      widk(1,cas)=max(data)-min(data);
    end;
    plus=(1./widk-1)*(n-1);
    plus=sort(plus);
    ptable(k,:)=plus(frac);
    wid(k,:)=widk;
    n,showtime
  end;
  n_plus=[nlist'/1000,ptable]
  return;
end;


data=data(:)'-min(data);
maxdata=max(data);
N=size(data,2);
if maxdata==0 | N<3,
  disp('data are constant or uninformative');
  teststat=NaN;
  testcor=NaN;
  return;
end;
data=sort(data)/maxdata; % normalized to range [0,1]



if 1,
  % just a few test statistics
  % found useful in this order
  teststat=zeros(1,3);
  u=[1:N];
  dev=(N-1)*data-u;
  teststat(1)=max(dev)-min(dev);
  e=round(0.75*N);
  teststat(2)=((N-1)/e)*max(data(1+e:end)-data(1:end-e));
  e=round(sqrt(2*N));
  dev=data(1+e:end)-data(1:end-e);
  teststat(3)=max(dev)-min(dev);
  dev=data(2:end)-data(1:end-1);
  teststat(4)=(N-1)*sum(dev.^2);
  testcor=[];
  return
end;


% old test

nt=40;
if N==3, ne==0;
else     ne=N-2;min(fix(N/2-1),50); 
         % just for beginning; later use powers of 2
end;
teststat=zeros(1,ne+1);
u=[1:N];
dev=(N-1)*data-u;
teststat(1)=max(dev)-min(dev);
for e=1:ne,
  dev=data(1+e:end)-data(1:end-e);
  teststat(e+1)=max(dev)-min(dev);
end;
testcor=[];
return


% old test

nt=40;
teststat=zeros(1,2*nt);
u=[1:N];
for i=1:nt,
  dev=(N-20+i)*data-u+(21-i)/2;
  teststat(i)=max(dev)-min(dev);
  teststat(i+nt)=sum(abs(dev));
end;
testcor=[];
return



% old test

teststat=zeros(1,6);
u=[1:N]/(N+1);
dev=data-u;
teststat(1)=max(dev);
teststat(2)=max(-dev);
teststat(3)=teststat(1)+teststat(2);
teststat(4)=max(dev.^2./(u.*(1-u)));

dif=(N+1)*(data(2:end)-data(1:end-1));
teststat(5)=max(dif);

curv=dif(2:end)-dif(1:end-1);
teststat(6)=max(curv);
teststat(7)=max(-curv);
teststat(8)=teststat(6)+teststat(7);

if N==3, ne==0;
else     ne=min(fix(N/2-1),15); 
         % just for beginning; later use powers of 2
end;
testcor=zeros(ne,5);
for e=1:ne,
  d1=dif(1:end-2*e);
  d2=dif(1+e:end-e);
  d3=dif(1+2*e:end);
  testcor(e,1)=max(dif(1+e:end).*dif(1:end-e));
  testcor(e,2)=max(d2.*(d1+d3));
  testcor(e,3)=max(d2.*max(d1,d3));
  testcor(e,4)=max(d2.*d1.*d3);
  testcor(e,5)=max(abs(curv(1+e:end).*curv(1:end-e)));
end;
  


