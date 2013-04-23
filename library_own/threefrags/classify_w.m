suffs = suffdata.suffs;
groupcount = suffdata.groupcount;
target = zeros(size(suffs));
nonempty = suffdata.nonemptygrps;

% initial classes
if 1,
    % most populated
    [n, perm] = sort(-groupcount);    
    ncl = 40;
    initcl = perm(1:ncl);
    target(initcl) = 1:ncl;    
    
else
    % random, as suggested
    % get random aa combinations in 20x20
    ncl = 40;
    p = randperm(400);   
    target(p(1:ncl)) = 1:ncl;
end

options.max_iterations = inf;
options.plots = 0;
options.ask = 0;
options.term_percentage = 0.01;
options.entweight = 0;
options.sort = 1;
options.merge = 1;

dprintf('maximum of %i groups, %i nonempty, %i initial classes', length(suffs), length(nonempty), max(target));

[potc,cl,confus,sufffreq,datafreq] = suffcov(suffs(nonempty),target(nonempty),options);

classanal;