function [zo,alphao,a] = makezforopticov(geom,all_s,alphas)
% I think z = tan((phi_i - alpha_i)/2), phi_i = tors(:,3) * pi/100, t_i =
% tors(:,3)
   % format of alphas:
   % ci cpi t1alpha t2alpha t3alpha (inf means only two classes)
% all_s = HEL27_anal08(geom);  list of s values for all 4-fragments, I
% think
% the name zo indicates that it is the original z, not the z * 100/pi
% (which is never used, anyway)
   
alphaz = alphas(:,3:5)*pi/100;
      
% don't seem to have to worry about the inf's - those combinations 
% (i.e. c = 1, cp = 1, t = 3, or c = 2, cp = 2, t = 3,... don't occur in a)
      
phio = double(geom(:,3))*pi/100;
a = sys10toX(all_s,[3 3 3]);
alphao = zeros(size(all_s));
zo = zeros(size(all_s));

% z depends on alpha_i
for i = 1:size(all_s,1)
  alphao(i) = alphaz(a(i,2) + (a(i,1) - 1)*3,a(i,3)); 
  zo(i) = tan(((phio(i)-alphao(i))/2));
end

