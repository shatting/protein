function [data,class] = suff_numrecon( suff, dataknown, idxknown )
%SUFF_NUMRECON Reconstruction of partial data from a covariance model
%given by sufficient statistic suff. 
%
% See also: numrecon.m

pot=suff_num2pot(suff, 0); ncl=length(pot); n=size(dataknown,1);
d=length(pot(1).mean); data=zeros(n,d,ncl); potval=zeros(n,ncl);

idxrecon = setdiff(1:d,idxknown);

for i=1:ncl
    
    mu=pot(i).mean; C=pot(i).cov;   
    data(:,:,i) = numrecon(mu',C,dataknown,idxknown);   
    y=data(:,:,i)'; MU=repmat(mu,1,n);
    potval(:,i)=dot((y-MU),inv(C)*(y-MU),1);

end

[void,class]=min(potval,[],2);

% vectorized selection of minimal value
sel=accumarray([(1:n)',ones(n,1),class],1,[n,1,ncl]);
data=sum(data.*repmat(sel,[1,d,1]),3);