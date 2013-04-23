function [ crossmat ] = crossmat( v )
%CROSSMAT a x b = crossmat(a)*b

error('use cross() instead of geom.crossmat');

crossmat(1,:) = [0, -v(3), v(2)];
crossmat(2,:) = [v(3), 0, -v(1)];
crossmat(3,:) = [-v(2), v(1), 0];