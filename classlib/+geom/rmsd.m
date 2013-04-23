function [ rmsd ] = rmsd( coords1,coords2 )
%RMSD Summary of this function goes here
%   Detailed explanation goes here

rmsd = sqrt(1/size(coords1,1)*sum(sum((coords1-coords2).^2)));

end

