function [ fval, individualpots, individualents ] = seq_pots(x, potseq, potential)
% [ fval, individualpots ] = seq_pots(x, potseq, potential)
% used for generating objective function as well as for evaluation of the elliptic
% constraints
%
% INPUT:
%       x           [c_1,..,c_(n-2),z_1,..,z_(n-3)]'
%       potseq      fragment index into potential struct array
%       potential   struct array (see opticov.m)
% OUTPUT:
%       fval        sum_i(V_{gamma,s}(x_{i:i+2} + ent_{gamma,s})
%       individualpots 

nfrags = length(potseq);

individualpots = zeros(nfrags,1)';
individualents = zeros(nfrags,1)';

for i = 1:nfrags,
    ind_c = 2*(i-1)+1 + [0 1];
    ind_z = nfrags + 1 + i;
    
    x_i = x([ind_c,ind_z]);
    pot_i = potential(potseq(i));
    
    Lxminusmu = pot_i.L*(x_i - pot_i.mean);
    individualpots(i) = Lxminusmu'*Lxminusmu;
    individualents(i) = pot_i.ent;
end

fval = sum(individualpots + individualents);

end
