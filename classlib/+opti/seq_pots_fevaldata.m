function [ fval, individualpots ] = seq_pots_fevaldata(obj,x)
% same as seq_pots/_fast
%  (0.938 secs->0.639 secs for all chains in db (fast))
%  (0.938 secs->0.134 secs for all chains in db (turbo))         
%
% INPUT:
%       x           [c_1,..,c_(n-2),z_1,..,z_(n-3)]'
%              
% OUTPUT:
%       fval        sum_i(V_{gamma,s}(x_{i:i+2} + ent_{gamma,s})
%       individualpots 

    xv = x(obj.idx);
    xminusmu = (xv - obj.means)';
    individualpots = obj.Ls.*xminusmu(obj.idx2,:);
    individualpots = sum(individualpots,2); % 
    individualpots = reshape(individualpots,3,obj.nfrags);
    individualpots = sum(individualpots.^2,1);

    fval = sum(individualpots + obj.individualents);

    if (nargout > 1 && (nargin < 3 || givegrad))
        %compute gradient

    end

end


