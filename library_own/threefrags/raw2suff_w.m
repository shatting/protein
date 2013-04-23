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

wtfn = @(x)max(0,(rawdata.rbins(end)^2-x.^2).^3);
maxwt = wtfn(rawdata.rbins(1));
wtfn_norm = @(x)(maxwt-wtfn(x))./maxwt;

wtfn_neum = @(x)wtfn(x)./(1+wtfn(x));

% groups formed by aa numbers of centers
grps = sysXto10(rawdata.istore(:,[2,5]), [20,20]);

nonempty = unique(grps);
groupcount = zeros(20*20,1);
suffs = cell(size(groupcount));
emptysuff = suffstat(1,zeros(size(feats,2),1));

dprintf('%i nonempty of %i groups',length(nonempty),20*20);

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
    r = feats(idx',end);
    
    wt=wtfn_norm(r);

    %wt2 = wtfn_neum(r);
    
    suffs{s} = suffstat(emptysuff,f,ngr,ones(ngr,1),[],wt);
    suffs{s}.group = s;    
end

suffdata.suffs = suffs;
suffdata.groupcount = groupcount;
suffdata.nonemptygrps = nonempty;
suffdata.raw_groups = grps;
suffdata.raw = rawdata;

figure;
hist(groupcount,100);
title 'group sizes';

