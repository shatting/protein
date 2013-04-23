function [ potval ] = pointpotential( pot, x )
%POINTPOTENTIAL x(:,i) is the i-th data vector

if (size(x,2) > 1)
    Lxminusmu = pot.L*(x - repmat(pot.mean,1,size(x,2))); % was non-vectorized: pot.L*(x - pot.mean);
    potval = sum(Lxminusmu.*Lxminusmu,1)'; % was non-vetorized: Lxminusmu'*Lxminusmu;
else
    Lxminusmu = pot.L*(x - pot.mean);
    potval = Lxminusmu'*Lxminusmu;
end

end

