

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% peakloc.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Npeak,mu,L]=peakloc(Npeak,data,mu,L);
% find peak given by ||(y-mu)*L||_2<=1 
% near approximate input peak ||(y-mu)*L||_2<=const
% note that data points are rows!
%
% Npeak           number of points in support of approximate peak
% data(l,:)	  l-th data point
% mu(1,:)	  approximate peak center 
% L        	  approximate scaling matrix (optional, L=1 allowed)
% 
% Npeak           number of points in support of peak
% mu(1,:)	  peak center 
% L        	  scaling matrix
%
function [Npeak,mu,L]=peakloc(Npeak,data,mu,L);

prt=0;    % printlevel (0=print nothing)

if nargin<4, L=1; end;

% transpose if necessary
d=length(mu);
if d<2, error('data must have at least 2 features'); end;
if size(mu,1)>1, mu=mu'; end;
if size(mu,1)>1, error('mu must be a vector'); end;
if size(data,2)~=d, data=data'; end;
if size(data,2)~=d, d, size_data=size(data), error('bad data size'); end;

N=size(data,1);
if N<2, error('too few data points'); end;
on=ones(N,1);
ind=[1:N];
it=0;
while 1,
  it=it+1;
  if prt, disp(['iteration ',num2str(it)]); end;
  e=round(sqrt(N));
  f=min(Npeak,fix(0.7*N));
  if ~(f>e & f<N-e), % needed for head and tail
    decrease=0; % small peak support or small sample; stop
  else 
    % compute potential
    Y=(data-mu(on,:))*L;
    V=0;
    dY=size(Y,2);
    for k=1:dY,
      V=V+Y(:,k).^2;
    end;
    % find peak in the good half of the V-distribution
    VV=sort(V);
    density=(2*e)./(VV(1+2*e:N)-VV(1:N-2*e));
    [head,peak]=max(density(1:f-e));
    tail=mean(density(f+1-e:end));
    thresh=peak+2*e+min(find(density(peak+1:end)<=tail));
    if prt,
      subplot(2,1,1);plot(VV);
      subplot(2,1,2);plot(density);
      Npeak,e,f-e,peak,thresh,input('cont>'); 
    end;
    Nold=Npeak;
    if isempty(thresh),
      decrease=0; 
      Vlim=VV(end);
    else
      Vlim=VV(thresh);
      ind=find(V<=Vlim);
      Npeak=length(ind);
      decrease=(Npeak<Nold);
    end;
  end;

  if prt, Npeak, end;
  if decrease | it==1,
    % compute potential for new peak support
    mu=mean(data(ind,:),1);
    Y=[data(ind,:)-mu(ones(length(ind),1),:)]';
    if norm(Y,inf)==0, L=1; return; end;
    [U,Sig]=svd(Y,0);
    % regularization
    sig=sqrt(diag(Sig));sig=sig/max(sig);
    D=diag(1./max(sig,0.05));
    condit=cond(D);
    L=U*D;
  else
    % smallest reasonable peak support reached
    L=L/sqrt(Vlim);
    return;
  end;

end;
  
  


