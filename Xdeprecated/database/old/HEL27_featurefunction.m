function [ features ] = HEL27_featurefunction(geomstruct,idx)
%HEL27_FEATURES Summary of this function goes here
%   Detailed explanation goes here

if (nargin == 0)
    features = 3;
    return;
else
    features = zeros(length(idx),3);
end
 
if (nargin == 1), % all fragments
    idx = 1:length(geomstruct.seq)-3;
end

features(:,1) = geomstruct.c(idx);
features(:,2) = geomstruct.c(idx+1);
features(:,3) = geomstruct.beta(idx,:);

end

