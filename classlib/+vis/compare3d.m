function [ ] = compare3d( obond, pbond ,opt)
% plots calpha coordinates of two proteins
% if input already in calphas, use opt = 1
% otherwise use opt = 0 or nothing


if nargin < 3
    opt = 0;
end

subplot(1,2,1);
if opt == 0
ribbon3d(bond2coords(obond),5);
else
ribbon3d(obond,5);
end

subplot(1,2,2);
if opt == 0
    ribbon3d(bond2coords(pbond),5);
else
    ribbon3d(pbond,5);
end
%also try mirrorz