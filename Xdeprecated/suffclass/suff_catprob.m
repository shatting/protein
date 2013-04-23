function [ pclfeats,pcl,pfeats,pfeatscl, direct ] = suff_catprob( suff, feats )
%SUFF_CATPROB Summary of this function goes here
%   Detailed explanation goes here

ncl = suff.ncl;
count=suff.count;

d = size(count,1);

Nclass = sum(suff.count,2);
Nclass = squeeze(Nclass(1,1,:));
Ndata = sum(Nclass);

for cl = 1:ncl, % calculate P(feats|cl)
    pfeatscl(cl,1) = prod(diag(count(1:d,feats,cl))/Nclass(cl));
end

% calculate P(feats)
pfeats = sum(pfeatscl);

% calculate P(cl)
pcl = Nclass/Ndata;

% calculate P( cl | feats )
pclfeats = pcl.*pfeatscl/pfeats;

%direct = (Ndata./Nclass).^(d-1) * 

end
