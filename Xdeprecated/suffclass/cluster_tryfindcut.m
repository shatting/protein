% CLUSTER_TRYFINDCUT Try to find a good cut in cluster merge sequence.
% [ cut ] = cluster_tryfindcut( cl, cost, options)

function [ cut ] = cluster_tryfindcut( cl, cost, options )

ncl = max(cl(:,end));
numclbds(1) = ceil(options.clust_minclfact*ncl);
numclbds(2) = ceil(options.clust_maxclfact*ncl);
relcostbds = options.clust_relcostbds;
maxdiff = options.clust_maxdiff;
minprev = options.clust_minprev;
minstep = options.clust_minstep;

cost = cost/max(cost);

diff = diffq(cost); % get difference quotients
diff = [-inf diff -inf]; % make same length as cost
nsplits = length(cost);
idx = 1:nsplits;

% xy bounds
numbds = nsplits - numclbds + 1;
numbds = numbds(2:-1:1);
xybds=cost>=relcostbds(1) & cost <= relcostbds(2) & idx>=numbds(1) & idx<=numbds(2);
diff(~xybds)=-inf;

% maxdiff bounds
diffbds = diff<=maxdiff;

steps = cost(2:end)-cost(1:end-1);
steps = [steps 0];
diffbds(steps<minstep) = 0;

% prev bounds
prev = find(cost(2:end)-cost(1:end-1) < minprev);
diffbds(prev+1) = 0;

% apply diffbds
diff(~diffbds) = -1;

[m,cut] = max(diff);

% save some info before possible error
cluster_info.numbds = numbds;
cluster_info.numclbds = numclbds;
cluster_info.diff = diff;
cluster_info.diffbds = diffbds;
 
if m<=0,
    warning('clustering not applied, parameters too rigid, please retry.');
else      
    
end

end