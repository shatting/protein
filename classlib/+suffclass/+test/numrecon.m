n = 100; % number of data points
d = 4;   % dimension
ndknown = 2; % number of known components
gauss = 0;   % normally distributed?
cestimate = 1; % estimate covariance?
gendata = 1;

if (gendata)
% generate data
    C = (rand(d)-0.5)*2;
    C = C*C'
    eigs = eig(C)
    mu = 10 * (rand(1,d)-0.5);
    data = gaussian(mu,C,n,gauss);

    % select known components
    idxknown = randperm(d);
    idxknown = sort(idxknown(1:ndknown));
    disp('known dimension(s):');
    disp(idxknown);
end

% reconstruct
if cestimate,
    suff = suffstat(1,zeros(d,1));
    suff = suffstat(suff,data,n,ones(n,1));
    datarecon = suff_numrecon(suff,data(:,idxknown),idxknown);
    scatterplot2(data,datarecon,idxknown, mu, C, suff);
else
    datarecon = numrecon(mu,C,data(:,idxknown),idxknown);
    scatterplot2(data,datarecon,idxknown, mu, C);
end
