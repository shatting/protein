

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dynoptR.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [p,y,Vmin]=dynoptR(T,V,U);
% finds global minimum of the potential encoded in dpot:
% min_{u,x} V(u,x):=sum_{t=1:n} V_t(u_t,x_t)+U(u_{n+1})
% s.t.      T_t(u_t,x_t)=u_{t+1}
% 
% T(u,x,t)        T_t(u,x), transition table (optionally only a matrix)
% V(u,x,t)        V_t(u,x), recurrent potential
% U(u,1)          U(u), optional terminal potential (default 0)
%
% p(1,:)          solution vector p
% y(1,:)          solution vector y 
% Vmin            optimal potential value
%
function [p,y,Vmin,U]=dynoptR(T,V,U)

[nu,nx,n]=size(V);
p=zeros(1,n+1);
y=zeros(1,n);
if nargin==2, Ut=zeros(nu,1); end;
Y=zeros(nu,n); % intermediate optimal points
U = zeros(nu,n); %mins from backwards loop, eventually need to take into account s_(n+1)

if ndims(T)==2,
  t=n;
  % backward loop
  while t>0, % reduce potential
    Vt = V(:,:,t)+ Ut(T);% repmat(Ut,1,9); 
    [U(:,t),Y(:,t)]=min(Vt,[],2); %finds smallest values in each row of Vt - U, indices of minimum values - X(:,t)
    Ut = U(:,t);
    t=t-1;
  end;
  %% check:
  %actmin = zeros(nu,1);
  %for i = 1:n
  %    actmin = actmin + (min(V(:,:,i)'))';
  %end
  %disp('actual sum of minimums should be equal to U(:,1) (if constraint is fulfilled)');
  %actmin == U(:,1)
  %actmin
  %actmin == W(:,end)
  %and the forward loop -- not actually needed for W, but just to see if
  %everything is working
  [Vmin,p(1)]=min(U(:,1)); %Vmin is min of U, so min of sum(V(:,:,t)), s(1) is index of Vmin in U -- starting value for forwards loop
  %[Vminf,r(n+1)] = min(W(:,n));
  % forward loop
  if size(T,1) > 1 %doesn't work for T = 1, but the min sequence is just Y in that case
      while t<n, % collect solution
        t=t+1;
        y(t)=Y(p(t),t); 
        p(t+1)=T(p(t),y(t));
      end;
  else
      y = Y;
  end
  
else %ndimsT > 3
     %I haven't changed this part yet, our T is 1-D
  t=n;
  % backward loop
  while t>0, % reduce potential
    Vt=V(:,:,t)+U(T(:,:,t));
    [U,X(:,t)]=min(Vt,[],2); %finds smallest values in each row of Vt - U, indices of minimum values - X(:,t)
    t=t-1;
  end;
  [Vmin,h(1)]=min(U);
  % forward loop
  while t<n, % collect solution
    t=t+1;
    x(t)=X(h(t),t);
    h(t+1)=T(h(t),x(t),t);
  end;
end;
