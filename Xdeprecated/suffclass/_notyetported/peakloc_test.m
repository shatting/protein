


ncas=1000;
Nlist=[100,300,1000,3000,10000,30000];
dlist=[500,200,100,50,20,10,5,3,2];

Nlist=[3000];
dlist=[2 3 5 10 20 50 100 200 500]
ncas=200;


NN=[];
if 1, 
  tic
  next=0;
  clear rr dd nn
  NN=Nlist;
end;

for N=NN,
  for d=dlist,
    next=next+1;
    res=zeros(1,ncas);
    for cas=1:ncas,
      if rem(cas,10)==0, disp(cas), end;
      % data=rand(N,d)-0.5;
      % data=[randn(N,d);randn(N,d)+sqrt(12/d)]; % distance sqrt(12)
      data=randn(N,d); % single Gaussian
      plot(data(:,1),data(:,2),'+','markersize',1)
      Npeak=N;
      mu=zeros(1,d);
      % mu=rand(1,d);
      [res(cas),mu,L]=peakloc(Npeak,data,mu);
    end;
    res=sort(res);
    ind=round([0.50,0.20,0.10,0.05,0.02,0.01]*ncas);
    rr(next,:)=res(ind);
    dd(next,1)=d;
    nn(next,1)=N;
    clf;plot(res/N)
    d_N=[d,N],showtime,disp(' ');
    drawnow;
    quot=rr(:,4)./nn
    % input('next>');
  end;
end;


yy=rr(:,4);
quot=yy./nn
clf;
subplot(1,2,1);plot(nn,yy./nn,'o')
subplot(1,2,2);plot(dd,yy./nn,'o')

return;




% old stuff


nlist=rsort(nn);
for k=1:length(nlist),
  n=nlist(k)
  ind=find(nn==n);
  [a(k,1),b(k,1),c(k,1)]=ratfit(dd(ind),yy(ind));
  % yy = a+b./(dd+c)
  clf;plot(dd(ind),yy(ind),'o');
  hold on;plot(dd(ind),a(k)+b(k)./(dd(ind)+c(k)),'-');
  hold on;plot(dd(ind),a(k)+1.3./(dd(ind)-1.3),':');
  input('next>');
end;
[nlist/1000,a,b,c]

return;

clf;plot(nn,yy-1.4-1.2./(dd-1.3)+6./sqrt(nn),'o')
clf;plot(dd,yy+6./sqrt(nn),'o')
[a,b,c]=ratfit(dd,yy+6./sqrt(nn))



