function [ features ] = HEFG_featurefunction( geomstruct,idx )
%HEFG_FEATURES = (c,cp,t,t')

if (nargin == 0)
    features = 8;
    return;
else
    features = zeros(length(idx),8);
end

if (nargin == 1), % all fragments
    idx = 1:length(geomstruct.seq)-3;
end
    
features(:,1) = geomstruct.c(idx);
features(:,2) = geomstruct.c(idx+1);
features(:,3) = geomstruct.t(idx,:);
features(:,4) = geomstruct.tp(idx,:);

end

