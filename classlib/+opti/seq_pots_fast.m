function [ fval, individualpots, individualents ] = seq_pots_fast(x, nfrags, idx, idx2, means, Ls, individualents)
% [ fval, individualpots ] = seq_pots_fast(x, potseq)
% same as seq_pots, but roughly 7 times(turbo)/1.45 times(fast) as fast. 
%  (0.938 secs->0.639 secs for all chains in db (fast))
%  (0.938 secs->0.134 secs for all chains in db (turbo))
% 
%
% INPUT:
%       x           [c_1,..,c_(n-2),z_1,..,z_(n-3)]'
%              
% OUTPUT:
%       fval        sum_i(V_{gamma,s}(x_{i:i+2} + ent_{gamma,s})
%       individualpots 
% if (nargin>2)
%     savefh(x);
% end

if 0, % fast verion 0.634, 349241 frags/s
    individualpots = zeros(3,nfrags);
    xv = x(idx);
    xminusmu = (xv - means);
    for i=1:nfrags,
       individualpots(:,i) = Ls(:,:,i)*xminusmu(:,i);
    end
    individualpots = sum(individualpots.^2,1);
else  % turbo version 0.134 sec, 1656311 frags/s   
    xv = x(idx);
    xminusmu = (xv - means)';
    individualpots = Ls.*xminusmu(idx2,:);
    individualpots = sum(individualpots,2); % 
    individualpots = reshape(individualpots,3,nfrags);
    individualpots = sum(individualpots.^2,1);
end

fval = sum(individualpots + individualents);

end
