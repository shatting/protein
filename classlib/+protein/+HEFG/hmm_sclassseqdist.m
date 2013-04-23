function [ d ] = hmm_sclassseqdist( sseq1, sseq2, distm)
%HMM_SCLASSSEQDIST 
% distm(i,j) 

test = 0;
if test,
    d2 = 0;
    idx = find(sseq1 ~= sseq2);
    for i=1:length(idx),
        d2 = d2 + distm(sseq1(idx(i)),sseq2(idx(i)));
    end
end

d = distm(sseq1,:)';
idx = 0:size(distm,2):size(distm,2)*(length(sseq1)-1); % starting index - 1
idx = idx' + sseq2;
d = d(idx);

d=sum(d);

if (test && d ~= d2)
    dprintf('error in dist calc.');
end

end

