

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% speclines.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function speclines(x,nr,range,spectext,sparse)
% prints spectral lines representing a 1D vector x
% with range range and ordinate [nr,nr+0.8] 
% text in spectext (default: none) is appended after the spectrum
% if sparse (default sparse=0), a random point on each spectral line 
%    is drawn instead of the whole spectral line, to see the density 
%    when there are too many points 
%
function speclines(x,nr,range,spectext,sparse)


if nargin<5, sparse=0; end;
if nargin>2,
  h=text(1.01*range(2)-0.01*range(1),nr+0.4,spectext);
  set(h,'fontsize',6);
end;
hold on;

n=length(x);
plot([range(1),range(1)],[nr,nr+0.8]);
plot([range(1),range(2)],[nr,nr]);
plot([range(1),range(2)],[nr+0.8,nr+0.8]);
plot([range(2),range(2)],[nr,nr+0.8]);

if sparse,
  y=nr+0.8*rand(n,1);
  plot(x,y,'.')
else
  for i=1:n,
    plot([x(i),x(i)],[nr,nr+0.8])
  end;
end;

