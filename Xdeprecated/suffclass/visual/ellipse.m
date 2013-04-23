% ELLIPSE Draw an ellipse.
% function h=ellipse(mu,r,C,lsty);
% draws an ellipse with center mu and covariance matrix r^2*C
% 
% lsty      line style (default '-b')
% h         line handle
%           set(h,'linestyle',sty)    changes the line style
%           set(h,'linewidth',wid)    changes the line width
%           set(h,'visible','off')    makes line invisible
%
function  h=ellipse(mu,r,C,lsty,col,wid);

if nargin==3, lsty='-b'; end;
if nargin<=5, wid=1; end;

R=chol(C);
phi=pi*[-50:50]/50;
x=r*R'*[cos(phi);sin(phi)];

if nargin<=4,
    h=plot(mu(1)+x(1,:),mu(2)+x(2,:),lsty);
else
    h=plot(mu(1)+x(1,:),mu(2)+x(2,:),lsty,'LineWidth',wid,'Color',col);
end
