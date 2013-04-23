function [z,alpharad] = HEL27_transform2z(geom,all_s)
% z = tan((beta_i - alpha_s{i})/2)
  
% get the alphas as they were used, despite the probable error
[bds, alphas] = HEL27_definition;
% format of alphas:
% 100*[ci cpi [t1alpha t2alpha t3alpha]/pi] (inf means only two classes)
   
% the following can be vectorized, somehow.
% but, as d. knuth put it: "premature optimization is the source of all evil"
alphasrad = alphas(:,3:5)*pi/100;
ccpt = sys10toX(all_s,[3 3 3]);
alpharad = zeros(size(all_s));
z = zeros(size(all_s));

for i = 1:length(all_s)
  alpharad(i) = alphasrad(ccpt(i,2) + (ccpt(i,1) - 1)*3, ccpt(i,3)); 
  z(i) = tan( (geom(i,3)-alpharad(i))/2 );
end

