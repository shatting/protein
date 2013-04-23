
data=randn(100,2,2); % 
%data(data>0.5)=data(data>0.5)+0.5;

u=min(data(:)); v=max(data(:)); d=1.3*(v-u);

clf;
for j=1:2,
 for k=1:2,
  x=data(:,j,k)+(k-1)*d;
  nr=-j;
  range=[u,v]+(k-1)*d;
  spectext=num2str(k);
  sparse=0;
  speclines(x,nr,range,spectext,sparse)
% prints spectral lines representing a 1D vector x
% with range range and ordinate [nr,nr+0.8] 
% text in spectext (default: none) is appended after the spectrum
% if sparse (default sparse=0), a random point on each spectral line 
%    is drawn instead of the whole spectral line, to see the density 
%    when there are too many points 
end;
end;
