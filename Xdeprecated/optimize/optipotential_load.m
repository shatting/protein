function [ optipot  ] = optipotential_load( probname )
%OPTIPOTENTIAL_LOAD Summary of this function goes here
%   Detailed explanation goes here

global optdatadir;

load([optdatadir,filesep,'optipot_',probname,'.mat']);



end

