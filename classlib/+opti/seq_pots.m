function [ fval, individualpots, individualents ] = seq_pots(x, potseq, savefh)
% [ fval, individualpots ] = seq_pots(x, potseq)
% used for generating objective function as well as for evaluation of the elliptic
% constraints
%
% INPUT:
%       x           [c_1,..,c_(n-2),z_1,..,z_(n-3)]'
%       potseq      struct array containing
% OUTPUT:
%       fval        sum_i(V_{gamma,s}(x_{i:i+2} + ent_{gamma,s})
%       individualpots 

if (nargin>2)
    savefh(x);
end

nfrags = length(potseq);

individualpots = zeros(nfrags,1)';
% want: [c_1,c_2,z_1 ; c_2,c_3,z_2 ; ..]'
idx = [(1:nfrags)',(2:nfrags+1)',(nfrags+2:2*nfrags+1)']';
xv = x(idx);
individualents = [potseq.ent];
xminusmu = xv - [potseq.mean];

%how to do next w/o loop?
for i=1:nfrags,
   Lxminusmu = potseq(i).L*xminusmu(:,i);
   individualpots(i) = Lxminusmu'*Lxminusmu;   
end
    
fval = sum(individualpots + individualents);

end
