function [ pot ] = opticov( feats, s, gamma )
%[ pot ] = opticov( feats, s, gamma )
% generates potential parameters for each (s,gamma) class combination

ns = max(s);
ng = max(gamma);

[n nf] = size(feats);

suff = suffstat(ns*ng, zeros(nf,1));

suff = suffstat(suff, double(feats), n, sysXto10([s gamma],[ns ng]));

pot = suff_num2pot(suff,1); % remove regularization? does not really make sense here

% works also, but get completely different values
%[meanx,Rx,beta,Vmin,sigx]=mom2pot(suff);
