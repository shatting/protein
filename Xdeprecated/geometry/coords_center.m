function [ centered_coords ] = coords_center( coords )
%[ centered_coords ] = coords_center( coords )

centered_coords = coords - repmat(mean(coords),size(coords,1),1);