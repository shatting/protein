
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% merge_bench.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation for analyzing significance test 
% for merging leaves in a decision tree
%
% Qmax=2.4 is adequate for s>25 to 100% and always to >96%
% Qmax=1+c/sqrt(s+1) is adequate for 
%  c    =  3  4  5  6    7
%  perc = 95 98 99 99.5 99.7


% takes 50 minutes on hektor
nall=[3000 1000 300 100 30 10];
sall=[2:4 6:2:20 25:5:50 60:10:100];
ncas=1000;

tic;
nn=length(nall);
nn2=nn*(nn+1)/2;
ns=length(sall);
n1list=zeros(ncas*ns*nn2,1);
n2list=zeros(ncas*ns*nn2,1);
slist=zeros(ncas*ns*nn2,1);
Qlist=zeros(ncas*ns*nn2,1);
table=zeros(ns*nn2,5);
now=0;
new=0;
for i1=1:nn,
  n1=nall(i1);
  for i2=1:i1,
    n2=nall(i2);
    N=n1+n2;
    clf;drawnow;hold on
    for s=sall,s
      new=new+1;
      f1=zeros(s,1);
      f2=zeros(s,1);
      QQ=zeros(ncas,1);
      for cas=1:ncas,
        if rem(cas,inf)==0, disp(cas); end;
        now=now+1;
        % random probability for all classes   
        sample1=rand(1,n1); 
        sample2=rand(1,n2); 
        bounds=zeros(s+1,1);
        for cl=1:s,
          p=rand;bounds(cl+1)=bounds(cl)+p;
        end;
        bounds=bounds/bounds(s+1);
        bounds(1)=-1;
        for cl=1:s, 
          f1(cl)=sum(sample1>bounds(cl)&sample1<=bounds(cl+1));
          f2(cl)=sum(sample2>bounds(cl)&sample2<=bounds(cl+1));
        end;
        f1=f1/n1;
        f2=f2/n2;
        S1=sum(f1.^2);
        S2=sum(f2.^2);
        S12=sum((f1-f2).^2);
        Q=[(n1-1)*(n2-1)*S12+(n1-1)*S2+(n2-1)*S1]/(n1+n2-2);
        QQ(cas)=Q;
        n1list(now)=n1;
        n2list(now)=n2;
        slist(now)=s;
        Qlist(now)=Q;
      end;
      Qmean=mean(QQ);Qstd=std(QQ);
      table(new,:)=[n1,n2,s,Qmean,Qstd];
      % plot(s+0*QQ,QQ,'o','markersize',1);
      % n1,n2,s
      % drawnow;
      showtime
    end;
    n1,n2,showtime
    disp('==========================================================');
    %input('next>');
  end;
end;
clf
subplot(2,2,1);
plot(n1list,Qlist,'o','markersize',1);set(gca,'xlim',[0,1100])
subplot(2,2,2);
plot(n2list,Qlist,'o','markersize',1);set(gca,'xlim',[0,1100])
subplot(2,2,3);
plot(n1list+n2list,Qlist,'o','markersize',1);set(gca,'xlim',[0,1100])

Qtable=[];
perc=[];
ii=0;
for s=sall,
  s
  ii=ii+1;
  ind=(slist==s);
  Q=sort(Qlist(ind));Q=Q(:)';
  nQ=length(Q);
  frac=round([0.9,0.95,0.99,0.999,1]*nQ);
  perc(ii,1)=max(find(Q<=2.4))*(100/nQ);
  perc(ii,2)=max(find(Q<=1+3/sqrt(s+1)))*(100/nQ);
  perc(ii,3)=max(find(Q<=1+4/sqrt(s+1)))*(100/nQ);
  perc(ii,4)=max(find(Q<=1+5/sqrt(s+1)))*(100/nQ);
  perc(ii,5)=max(find(Q<=1+6/sqrt(s+1)))*(100/nQ);
  perc(ii,6)=max(find(Q<=1+7/sqrt(s+1)))*(100/nQ);
  Qtable(ii,:)=Q(frac);
  % plot(Q);set(gca,'ylim',[0,5]);
  % input('next>');
end;
s_Q=[sall',Qtable]
subplot(2,2,4);
cla;
plot(slist,Qlist,'o','markersize',1);hold on;
plot(slist,1+5./sqrt(slist+1))
% plot(slist,2+3./(slist1))
plot(sall,Qtable,'+');hold off;
s_perc=[sall',perc]

return;

clf;
plot(1./sqrt(sall),Qtable,'+');
plot(sqrt(sall),Qtable,'+');

