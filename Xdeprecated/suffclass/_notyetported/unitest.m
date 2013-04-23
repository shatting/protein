

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% unitest.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [fac,perc,kuiper]=unitest(data);
% Kuiper test statistics for testing uniformity of data 
%
% roughly perc percent of uniform samples exceed the test statistic
% fac   1.7  1.8  1.9  2.0  2.1  2.2   2.3   2.4   2.5   2.6   2.7   
% perc  6.8  4.6  3.1  1.8  1.1  0.56  0.28  0.11  0.06  0.03  0.015
%
% thus fac<=1.8, giving
%   kuiper<1.8*(sqrt(N-1.5)-0.7),
% is a good uniformity test
%
function [fac,perc,kuiper]=unitest(data);

data=sort(data);
range=data(end)-data(1);
if range<=0,
  kuiper=0;fac=0;perc=100;
  return;
end;
N=length(data);
data=data*((N-1)/range)-[1:N]; 
kuiper=max(data)-min(data);
fac=kuiper/(sqrt(N-1.5)-0.7);
if nargout==1, return; end;
percent=exp(4.8+fac.*(1.5-1.8*fac));
