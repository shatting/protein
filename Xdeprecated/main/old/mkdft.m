function [ dft ] = mkdft( data )
%MKDFT Summary of this function goes here
%   Detailed explanation goes here
% data(i,:) i-th vector
% dft(y1,y2,..) number of occurences of vector (y1,y2,..)

ymax = double(max(data,[],1));
dft = zeros(ymax);

for i=1:size(data,1),
    vec = double(data(i,:));
    %idx = vec*ymax' + 1;
    dft(vec(1),vec(2),vec(3),vec(4),vec(5)) = dft(vec(1),vec(2),vec(3),vec(4),vec(5))+1;
end
