function [ data ] = numrecon( mu, C, dataknown, idxknown )
%SUFF_NUMRECON Reconstruction of partial data from a covariance model
%given by sufficient statistic suff. 
% INPUT: 
%       mu          ... 1xd, covariance model mean
%       C           ... dxd, covariance model covariance matrix
%       dataknown   ... nxk, the k known components of data vectors, k<=d
%       idxknown    ... 1xk, index set of the k known components in the
%                            d-dimensional full vectors
% OUTPUT:
%       data ... reconstruced data set. 
%                data(:,idxknown) = dataknown
%                data(:,1:d/idxknown) = reconstructed data

d = length(mu);
[n,k] = size(dataknown);

idxrecon = setdiff(1:d,idxknown);

C11 = C(idxknown,idxknown);
C21 = C(idxrecon,idxknown);

meanknown = mu(idxknown);
meanrecon = mu(idxrecon);

data = zeros(n,d);

epsknown = dataknown - repmat(meanknown,n,1);

datarecon = epsknown * inv(C11)' * C21' + repmat(meanrecon,n,1);

data(:,idxknown) = dataknown;
data(:,idxrecon) = datarecon;

end
