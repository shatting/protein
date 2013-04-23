function [ s ] = seq2sclasses( bond )
%SEQ2SCLASSES quick & dirty, bond die aus RT*.mat
% zb. seq2sclasses(data{1}.seq,data{1}.bond)

geom = bond2geometry(bond);

s = HEL27_anal08([geom.c geom.cp geom.phi])';
