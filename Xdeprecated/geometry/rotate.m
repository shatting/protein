function [ x ] = rotate( x, Q )
%ROTATE rotates the vectors x(i,:) using Q

x = x*Q';