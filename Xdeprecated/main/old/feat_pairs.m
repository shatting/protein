function [ pairfeats ] = feat_pairs( feats, dontcheck )
%[ pairfeats ] = feat_pairs( feats )
% eg:
% fragfeats = [   seq(:,1) + 20*(seq(:,2)-1),...
%                 seq(:,1) + 20*(seq(:,3)-1),...
%                 seq(:,1) + 20*(seq(:,4)-1),...       
%                 seq(:,2) + 20*(seq(:,3)-1),...
%                 seq(:,2) + 20*(seq(:,4)-1),...
%                 seq(:,3) + 20*(seq(:,4)-1)];   

if (nargin<2), dontcheck = false; end

[n,d] = size(feats);

pairfeats = zeros(n,d*(d-1)/2);

sys = max(feats,[],1);

% find suitable data type
if (~dontcheck)
    m = max(feat_pairs(double(sys),true));
    if m > 255,
        feats = uint16(feats);
        pairfeats = uint16(pairfeats);
    end
    if m > 255*255,
        feats = unit32(feats);
        pairfeats = uint32(pairfeats);
    end
    sys = max(feats,[],1); % convert sys to same type as feats
end

% make pairfeatures
k = 1;
for i=1:d-1,
   
    for j=i+1:d,
       pairfeats(:,k) = feats(:,i) + sys(i)*(feats(:,j)-1);
       k=k+1;
    end
    
end    