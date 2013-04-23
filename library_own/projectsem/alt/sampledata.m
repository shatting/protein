function [ sample ] = sampledata( data, numaa )
%TESTDATA Summary of this function goes here
%   Detailed explanation goes here

nchains = min(numaa,length(data));
perm = randperm(length(data));

j = 1;

perm = sort(perm(1:nchains));

sample = cell(nchains, 1);

for i=perm,
   sample{j}.res = data{i}.res;
   sample{j}.name = data{i}.name;
   sample{j}.seq = data{i}.seq;
   sample{j}.bond = data{i}.bond;
   
   j = j + 1;
end