suffs = suffdata.suffs;
groupcount = suffdata.groupcount;
target = zeros(size(suffs));
nonempty = suffdata.nonemptygrps;

% initial classes
if 0,
    % most populated
    [n, perm] = sort(-groupcount);    
    ncl = 80;
    initcl = perm(1:ncl);
    target(initcl) = 1:ncl;    
    
else
    % random, as suggested
    % get random aa combinations in 20x20
    naacombis = 9;
    allow = 1:400;
    seedcl = [];
    for i=1:naacombis,
        g = ceil(length(allow)*rand);
        aa = sys10toX(allow(g),[20,20]);
        aa = repmat(aa,nrbins,1);
        seedcl = [seedcl; aa [1:nrbins]'];
        allow(g) = [];
    end
    
    seedcl = sysXto10(seedcl,[20 20 nrbins]);
    
    target(seedcl) = 1:length(seedcl);
end

options.max_iterations = inf;
options.plots = 0;
options.ask = 0;
options.term_percentage = 0.01;
options.entweight = 0;
options.sort = 1;
options.merge = 1;

dprintf('maximum of %i groups, %i nonempty, %i initial classes', length(suffs), length(nonempty), max(target));

[potc,cl,confus,freq] = suffcov(suffs(nonempty),target(nonempty),options);

classanal;