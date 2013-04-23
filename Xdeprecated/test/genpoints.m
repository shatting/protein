function [ p ] = genpoints( mu, C, n )
%GENPOINTS generate n points normally distributed around mu with covariance
%C

dim = length(mu);

if size(mu,1) > size(mu,2),
    mu = mu'
end;

x = randn(n,dim);

p = x * C + repmat(mu,n,1);

%clf;hold on; c = rand(2,2)-0.5; plot_data_n_suff(randn(1000,2)*(c*c'));
%ellipse([0,0],0.5,c*c',':'); hold off;