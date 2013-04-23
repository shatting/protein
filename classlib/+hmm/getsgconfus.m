function [ sgconfus ] = getsgconfus( dataset, sclassifier, gclassifier )
%GETSGCONFUS Summary of this function goes here
%   Detailed explanation goes here

sgconfus = getfreq([gclassifier.classify(dataset), sclassifier.classify(dataset)]); % we want gamma in the rows  

end

