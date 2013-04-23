% rawdata.readme
% 
% ans =
% 
% only threefrag pairs with rbins(1)<dist(i+1,j+1)<=rbins(end)     
% dstore(l,:) = [gram(gramidx) dist(i+1,j+1)]  double              
%               [   1:10             11     ]                      
% istore(l,:) = [seq(i:i+2) seq(j:j+2) rbin ]  uint8               
%               [    1:3        4:6      7  ]                      
% generated with getrawdata.m      



suffdata = struct;

disp('making sufficient statistics');
nrbins = length(rawdata.rbins)-1;
feats = rawdata.dstore(:,1:11);

% groups formed by aa numbers of centers and the r bin
grps = sysXto10([rawdata.istore(:,[2,5]) rawdata.istore(:,end)], [20,20,nrbins]);

nonempty = unique(grps);
groupcount = zeros(20*20*nrbins,1);
suffs = cell(size(groupcount));
emptysuff = suffstat(1,zeros(size(feats,2),1));

dprintf('%i nonempty of %i groups',length(nonempty),20*20*nrbins);

i = 0;
n=length(nonempty);
report = round(n/10);
for s = nonempty,
    i = i+1;
    if ~mod(i,report), dprintf('group %i of %i',i,n); end
    idx = grps==s;
    ngr = sum(idx);
    groupcount(s) = ngr;
    f = feats(idx',:);
    suffs{s} = suffstat(emptysuff,f,ngr,ones(ngr,1));
    suffs{s}.group = s;    
end

suffdata.suffs = suffs;
suffdata.groupcount = groupcount;
suffdata.nonemptygrps = nonempty;
suffdata.raw_groups = grps;
suffdata.raw = rawdata;

figure;
hist(groupcount);
title 'group sizes';

