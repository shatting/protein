function [ geom ] = prediction2geom( prediction, sseq, alphas )
%[ geom ] = prediction2geom( prediction )
% inverse of makezforopticov.m
% INPUT: 
%       prediction      [c_1,..,c_n-2,z_1,..z_n-3]
%       sseq            

alphaz = alphas(:,3:5)*pi/100;
           
phio = double(geom(:,3))*pi/100;
a = sys10toX(all_s,[3 3 3]);
alphao = zeros(size(all_s));
zo = zeros(size(all_s));

% z depends on alpha_i
for i = 1:size(all_s,1)
  alphao(i) = alphaz(a(i,2) + (a(i,1) - 1)*3,a(i,3)); 
  zo(i) = tan(((phio(i)-alphao(i))/2));
end


end
