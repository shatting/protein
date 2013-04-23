function [ d ] = kldiv( fp, fq )
%KLDIV kullback-leibler divergence D(p||q)
% [ d ] = kldiv( fp, fq )
% fp(:,i) is the i-th distribution to be tested, given by frequencies
% fq(:,1) is the reference distribution
%
% d(1,i) is D(fp(:,i)||fq)

reg = 0.001;

dimp = size(fp,1);
np = size(fp,2);

countp = sum(fp);
countpx = repmat(countp,dimp,1);
nq = sum(fq);

prp = fp./countpx;
prq = (fq+reg)/(nq+reg*dimp); % regularize

prqx = repmat(prq,1,np);

ratio = prp./prqx + eps;

logratio = log2(ratio);

weightedlogratio = prp.*logratio;

d = sum(weightedlogratio,1);

% conf = getfreq([cl1 cl2]);
% [n1 n2] = size(conf);
% f1 = sum(conf,2);
% f2 = sum(conf,1);
% n = sum(f1);
% 
% p1 = f1/n;
% f2x = repmat(f2,n1,1);
% %p2 = f2/n;
% 
% p1x = repmat(p1,1,n2);
% %p2x = repmat(p2,n1,1);
% p2x = conf./f2x;
% 
% rat = p1x./p2x;
% 
% lo = log2(rat);
% 
% pro = p1x.*lo;
% 
% d = sum(pro,1);

end









