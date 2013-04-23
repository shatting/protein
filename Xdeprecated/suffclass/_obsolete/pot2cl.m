

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pot2cl.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [pred,ppot]=pot2cl(meanx,Rx,beta,Vmin,data,logpr);
% predict classes for data from potential
%    V(x,t)=beta(t)*(x-meanx(:,t))'*Rx(:,:,t)*(x-meanx(:,t))+Vmin(t)
% created by mom2pot.m
%
% For prediction with prior class probabilities pr(t),
% minimize V(t|x)=V(x,t)-log(pr(t)),
% i.e., use Vmin-log(pr) in place of Vmin
%
% meanx(:,t)      mean for class t (rounded to integer)
% Rx(:,:,t)       upper triangular quadratic term for class t
% beta(1,t)       potential scaling factor
% Vmin(1,t)       potential minimum
% data(l,:)       l-th data vector, = x^T
%
% pred(l,1)       predicted class of l-th data vector
% ppot(l,1)       selected potential of l-th data vector
%                 
function [pred,ppot]=pot2cl(meanx,Rx,beta,Vmin,data)

if isempty(meanx),
  disp('data set is too small');
  error(Rx); % Rx contains the error message
end;

ncl=size(beta,2);
N=size(data,1);
pred=zeros(N,1);
ppot=zeros(N,1);
if N>=1000, t0=toc; end;
for l=1:N,
  if rem(l,1000)==0, 
    t=toc-t0;
    disp([num2str(l),'/',num2str(N),' in ',showtime(t),...
          ' of ',showtime(t*N/l)]); 
  end;
  Vl=inf;cl=0;
  x=data(l,:)';
  for t=1:ncl,
    dx=x-meanx(:,t);
    Vt=beta(t)*dx'*Rx(:,:,t)*dx+Vmin(t);
    if Vt<Vl, cl=t;Vl=Vt; end;
  end;
  pred(l)=cl;
  ppot(l)=Vl;
end;


