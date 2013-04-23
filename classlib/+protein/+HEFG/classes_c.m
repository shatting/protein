function [ cl ] = classes_c( frags )
%INITIALCLASSIFIER Summary of this function goes here
%   Detailed explanation goes here

H=(frags(:,1:2)<=20);
E=(frags(:,1:2)>50);
L=(frags(:,1:2)>20 & frags(:,1:2)<=50);
N=~(H|E|L);
cl=H+2*L+3*E+100*N;