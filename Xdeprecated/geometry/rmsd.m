function [ rmsd ] = rmsd( x, y )
%RMSD simply compute the rmsd of 2 sets of vectors. 
% to find optimal transform, use register.m

rmsd = sqrt(1/size(x,1)*sum(sum((x-y).^2)));
