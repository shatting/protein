function [ data ] = suff_numrecon( suff, dataknown, idxknown )
%SUFF_NUMRECON Reconstruction of partial data from a covariance model
%given by sufficient statistic suff. 
%
% See also: numrecon.m

pot = suff_num2pot(suff, 0);

data = numrecon(pot.mean',pot.cov,dataknown,idxknown);

end
