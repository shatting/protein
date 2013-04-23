%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% unitest_bench.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% benchmarking for unitest
% it seems that no test is better than the 
% Kuiper variant of the K/S test
%



if 0,
  clear;clear mex
  tic
  % calibrate Kuiper's test
  ncas=20000;
  nrep=5;
  fac=[10:30]/10;
  nf=length(fac);
  Nlist=[2:5 10 20 50 100 200 500 1000 2000 5000];
  % Nlist=[5 10 20 50 100 200 500 1000]; % N<5 is slightly off
  NN=length(Nlist);
  quant=zeros(7,nrep,NN);
  perc=zeros(nf,nrep,NN);
  cas=1;
  for k=1:NN,
    N=Nlist(k);
    for rep=1:nrep,
      unif=zeros(ncas,1);
      for cas=1:ncas,
        [ufac(cas),uperc(cas),unif(cas)]= unitest(rand(N,1));
      end;
      Q=sort(unif);
      nQ=length(Q);
      frac=round([0.5,0.7,0.9,0.95,0.99,0.999,1]*nQ);
      quant(:,rep,k)=Q(frac);
      Qf=sort(ufac);
      nQ=length(Qf);
      for j=1:nf,
        perc(j,rep,k)=100-max(find(Qf<=fac(j)))*(100/nQ);
      end;
    end;
    N,showtime
  end;

  fac=fac(1:end-2);
  pp=max(max(perc(:,:,4:NN),[],2),[],3);
  lp=log(pp(1:end-2))';
  [coef,tol]=polfit(fac,lp,0,2)
  clf;
  plot(fac,lp,'+');hold on
  plot(fac,4.8+fac.*(1.5-1.8*fac))
  for j=1:length(fac),
    pe=perc(j,:,1:NN);pe=pe(:);pe=pe(pe>0);
    plot(fac(j),log(pe),'.','markersize',3)
  end;
  hold off
  return

end;





if 0,
  % test bump-dip family
  for alp=0.5:9.5,
    sub=0;
    for lam=-0.5:0.1:0.51,
      lam=round(10*lam)/10;
      sub=sub+1;
      N=1000;
      x=alp*(sort(rand(N,1))-0.5);
      y=x-lam;;
      ind=(x>=1);y(ind)=x(ind)+lam;
      ind=(abs(x)<1);xind=x(ind);
      y(ind)=xind.*(1+lam*(2+xind.^2.*(abs(xind)-2)));
      y=y/max(abs(y));
      subplot(3,4,sub);
      % plot(x,y);
      hist(y,10);set(gca,'xlim',[-1,1]);
      xlabel(['a=',num2str(alp),' l=',num2str(lam)]);
    end;
    input('next>');
  end;
  return;
end;


if 1,
  % test a few statistics only
  ncas=100;
  alplist=[ 0.5  0.5  0.5  0.5  1    1    1    1  8    8    8    8];
  lamlist=[-0.5 -0.25 0.25 2   -0.5 -0.25 0.25 2 -0.5 -0.25 0.25 2];
  lamlist=1;
  alplist=2;
  symlist='x+*osdv^<>ph'; % 12 classes
  nc=length(alplist);
  nl=length(lamlist);
  if nl~=nc,error('sizes do not match'); end;
  N=100;
  for cas=1:ncas,
    if rem(cas,10)==0, disp(cas); end;
    % uniform data
    data=rand(N,1);
    teststat=unitest_try(data);
    test=teststat(:)';
    if cas==1, 
      nf=length(test); 
      unif=zeros(ncas,nf); 
      nonunif=zeros(ncas,nf,:);
    end;
    unif(cas,:)=test;

    % alternative data from various classes
    for cl=1:nc,
      alp=alplist(cl);
      lam=lamlist(cl);
      x=alp*(sort(rand(N,1))-0.5);
      y=x-lam;;
      ind=(x>=1);y(ind)=x(ind)+lam;
      ind=(abs(x)<1);xind=x(ind);
      y(ind)=xind.*(1+lam*(2+xind.^2.*(abs(xind)-2)));
      data=y/max(abs(y));
      teststat=unitest_try(data);
      test=teststat(:)';
      nonunif(cas,:,cl)=test;
    end;
  end;
  
  ms=5;
  for cas=1:3,
    figure(cas);clf;
    if cas==1,     j=1;k=2;
    elseif cas==2, j=1;k=3;
    else           j=2;k=3;
    end;
    clf;hold on
    for cl=1:nc,
      sym=symlist(cl);
      plot(nonunif(:,j,cl),nonunif(:,k,cl),sym,'markersize',ms);
    end;
    plot(unif(:,j),unif(:,k),'ro','markersize',5);hold off
  end;
  return;
end;




% look for fixed N on unitest_try.n for alternatives
ncas=300;
lamlist=[-0.5  -0.2  0.2  0.5  1  5];
alplist=[0.5  1  2  4  8];
%lamlist=[-0.3  -0.1  0.1  0.3  1  3];
%alplist=[0.5  1  2  4  8];
nl=length(lamlist);
na=length(alplist);
nc=nl*na;
for N=200,
  for cas=1:ncas,
    if rem(cas,10)==0, disp(cas); end;
    % uniform data
    data=rand(N,1);
    [teststat,testcor]=unitest_try(data);
    test=[testcor(:);teststat(:)]';
    if cas==1, 
      nf=length(test); 
      unif=zeros(ncas,nf); 
      nonunif=zeros(ncas,nf,nc);
    end;
    unif(cas,:)=test;

    % alternative data from various classes
    cl=0;
    for lam=lamlist,
      for alp=alplist,
        cl=cl+1;
        x=alp*(sort(rand(N,1))-0.5);
        y=x-lam;;
        ind=(x>=1);y(ind)=x(ind)+lam;
        ind=(abs(x)<1);xind=x(ind);
        y(ind)=xind.*(1+lam*(2+xind.^2.*(abs(xind)-2)));
        data=y/max(abs(y));
        [teststat,testcor]=unitest_try(data);
        test=[testcor(:);teststat(:)]';
        nonunif(cas,:,cl)=test;
      end;
    end;
  end;
end;



%clf
nrow=11;ncol=10;
kbest=1;
frac=0.95; 
unisort=sort(unif(:,kbest));
thresh=unisort(round(frac*ncas));
times=zeros(nf,1);
th=zeros(nf,1);
for k=1:nf,
  %subplot(nrow,ncol,k)
  unisort=sort(unif(:,k));
  xlim=[unisort(1),unisort(round(frac*ncas))];
  th(k)=xlim(2);
  %plot(unisort,zeros(ncas,1),'x','markersize',2);hold on
  tim=0;
  for cl=1:nc,
    ind=find(nonunif(:,kbest,cl)<thresh);
    tim=tim+sum(nonunif(ind,k,cl)<=th(k));
    %plot(nonunif(ind,k,cl),cl+0*ind,'x','markersize',2);
  end;
  times(k)=tim;
  %if ~any(isnan(xlim)), set(gca,'xlim',xlim,'xtick',[]); end;
  %set(gca,'ylim',[-1,nc+1],'ytick',[]);
  %xlabel([num2str(k)]);
  if rem(k,ncol)==0, drawnow; end;
  % input('next>');
end;


k_tim=[[1:nf]'  times/times(kbest)-1]
times(kbest)

plot(times/times(kbest)-1,'o'); hold on






return;



