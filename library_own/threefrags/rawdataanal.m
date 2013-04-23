
if (~exist('rawdata','var'))
    load rawdata;
end

ndata = size(rawdata.dstore,1);

% get things
rnorms = rawdata.dstore(:,11);          %dist(i+1,j+1)
rbins = rawdata.istore(:,end);


% interesting: has peak at ~6.3, then falls off again
subplot(2,1,1);
hist(rnorms,20);
title 'original r distances';
subplot(2,1,2);
hist(double(rbins));
title 'r bins';

