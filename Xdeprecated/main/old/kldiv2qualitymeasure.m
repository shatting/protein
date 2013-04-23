function [ wmeankld, meankld ] = kldiv2qualitymeasure( target, cl, varargin )
%KLDIV2QUALITYMEASURE Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(varargin)
    method = varargin{1};
else
    method = [];
end

conf = getfreq([target cl]);

p = conf./repmat(sum(conf,1),size(conf,1),1) + eps;
q = sum(conf,2)/sum(sum(conf)) + eps;

% test ps against qs
for i=1:size(p,2),
    m(i) = kldiv2(p(:,i),q,method);
end

p = sum(conf,1)/sum(conf(:));

meankld = mean(m);
wmeankld = m*p';

end
