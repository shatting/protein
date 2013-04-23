%GENPOINTS Generate points normally distributed around a mean with
%covariance.
function p = gaussian( mu, C, n, normal )

dim = length(mu);

if size(mu,1) > size(mu,2),
    mu = mu'
end;

if nargin < 4 || normal
    x = randn(n,dim);   
else
    x = rand(n,dim)-0.5;
end
p = x * C' + repmat(mu,n,1);

%clf;hold on; c = rand(2,2)-0.5; plot_data_n_suff(randn(1000,2)*(c*c'));
%ellipse([0,0],0.5,c*c',':'); hold off;