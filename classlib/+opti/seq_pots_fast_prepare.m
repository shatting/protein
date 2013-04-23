function [ nfrags, idx, idx2, means, Ls, individualents ] = seq_pots_fast_prepare( potseq )
%SEQ_POTS_INTERFACE Summary of this function goes here
%   Detailed explanation goes here

nfrags = length(potseq);
% want: [c_1,c_2,z_1 ; c_2,c_3,z_2 ; ..]'
idx = [(1:nfrags)',(2:nfrags+1)',(nfrags+2:2*nfrags+1)']';
individualents = [potseq.ent];
means = [potseq.mean];

if 0, % fast option
    Ls = zeros(3,3,nfrags);
    for i=1:nfrags,
        Ls(:,:,i) = potseq(i).L;
    end
    idx2 = [];
else % turbo option
    idx2 = 1:nfrags;
    idx2 = idx2([idx2;idx2;idx2]);    
    idx2 = reshape(idx2,3*nfrags,1);
    
    Ls = zeros(3*nfrags,3);
    for i=1:nfrags,
        Ls(3*(i-1)+1:3*i,:) = potseq(i).L;
    end
end

end

