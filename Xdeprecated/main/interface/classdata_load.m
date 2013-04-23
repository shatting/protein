function [ new_deal_data ] = classdata_load( problemname )
%CLASSDATA_LOAD Summary of this function goes here
%   Detailed explanation goes here
global classdatadir;
load([classdatadir,filesep,problemname]);

end

