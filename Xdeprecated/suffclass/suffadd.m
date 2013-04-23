% SUFFADD Add two sufficient statistics.
%
% function suff=suffadd(suff1,suff2,cas);
% add sufficient statistics suff1, suff2 and assign it to suff
% cas=1:     same classes
% cas=2:     create new class for each suff2 class

function s1=suffadd(s1,s2,cas)
 if nargin < 3 || cas==1,
   s1.moment=s1.moment+s2.moment;
   s1.count=s1.count+s2.count;
 else
   s1.ncl=s1.ncl+s2.ncl;
   s1.moment(:,s1.ncl+1:s2.ncl)=s2.moment;
   s1.count(:,s1.ncl+1:s2.ncl)=s2.count;
 end;
   
end