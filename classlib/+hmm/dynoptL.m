%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dynoptL.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [p,y,Vmin]=dynoptL(S,V,W);
% finds global minimum of the potential encoded in dpot:
% min_{u,x} V(u,x):=sum_{t=1:n} V_t(u_t,x_t)+W(u_{n+1})
% s.t.      S_t(u_{t+1},x_{t+1})=u_{t}
% 
% S(u,x,t)        S_t(u,x), transition table (optionally only a matrix)
% V(u,x,t)        V_t(u,x), recurrent potential
% U(u,1)          U(u), optional terminal potential (default 0)
%
% p(1,:)          solution vector p
% y(1,:)          solution vector y
% Vmin            optimal potential value
%
function [p,y,Vmin,W]=dynoptL(S,V,W);

[nu,nx,n]=size(V);
p=zeros(1,n+1);
y=zeros(1,n);
if nargin==2, Wt=zeros(nu,1); end;
Y=zeros(nu,n); % intermediate optimal points
W = zeros(nu,n); %mins from forwards loop
Wt = zeros(nu,1);

if ndims(S)==2,
  t=0;
  % forward loop
  while t<n, % reduce potential
    t = t + 1;
    Vt=V(:,:,t)+ Wt(S); %repmat(Wt,1,9); %
    [W(:,t),Y(:,t)] = min(Vt,[],2); %finds smallest values in each row of Vt - U, indices of minimum values - X(:,t)
    Wt = W(:,t); 
  end;
  
  % check:
  %actmin = zeros(nu,1);
  %for i = 1:n
  %    actmin = actmin + (min(V(:,:,i)'))';
  %end
  %disp('actual sum of minimums should be equal to W(:,end) (if constraint is fulfilled)');
  %actmin == W(:,end)
  [Vmin,p(n+1)]=min(W(:,n)); %Vmin is min of U, s(1) is index of Vmin in U -- starting value for forwards loop
  % backward loop
  % t = n
    if size(S,1) > 1 %doesn't work for T = 1, but the min sequence is just Y in that case
          while t>0, % collect solution
            y(t) = Y(p(t+1),t);
            p(t) = S(p(t+1),y(t));
            t = t-1;
          end;
    else
        y = Y;
    end
else %ndimsT > 3
    %I haven't changed this part yet, our T is 1-D
  t=n;
  % backward loop
  while t>0, % reduce potential
    Vt=V(:,:,t)+U(S(:,:,t));
    [U,X(:,t)]=min(Vt,[],2); %finds smallest values in each row of Vt - U, indices of minimum values - X(:,t)
    t=t-1;
  end;
  [Vminf,s(1)]=min(U);
  % forward loop
  while t<n, % collect solution
    t=t+1;
    x(t)=X(s(t),t);
    s(t+1)=T(s(t),x(t),t);
  end;
end;
