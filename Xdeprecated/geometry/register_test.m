function [ e, d] = register_test( n )
%RMSD_TEST Summary of this function goes here
%   Detailed explanation goes here

x = rand(n,3)*rand*100; %stretch
x = shift(x,rand(1,3)*rand*50); %translate
y = rand(n,3)*rand*100;
y = shift(y,rand(1,3)*rand*50);

[Q,q,rmr] = register(x,y);

rmd = rmsd(x,shift(rotate(y,Q),q));

e = (rmd-rmr)/rmd;
d = det(Q);