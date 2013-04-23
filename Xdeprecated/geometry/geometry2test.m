% round trip test

[x, alphas] = HEL27_definition;

% get a bond vector
maxb = -inf;
for p = 1:1334
    if mod(p,100) == 1
        p
        maxb
    end

bond = double(data{p}.bond);
bond = rotatematrix(bond); %bond./repmat(sqrt(sum(bond.^2,2)),1,3);

% get t geometry
geom = bond2geometry(bond);

% convert to z geometry
geomz = geometry2geometryz(geom, alphas);

% convert back to t geometry
geom_after = geometryz2geometry(geomz);

% convert back to bond
bond_after = geometry2bondn(geom_after);


bondn = bond./repmat(sqrt(sum(bond.^2,2)),1,3);

maxbond = max(max((bond_after - bondn)));
if maxbond > maxb
        maxb = maxbond;
    end
if maxbond > 10^(-5)
        problems = [problems,p];
end


% compare both geometries
% cmpstructs(geom,geom_after);
% 
% % compare both bonds
% % normalize original bonds first
%bondn = bond./repmat(sqrt(sum(bond.^2,2)),1,3);
% 
%[Q,q,rmse_bond] = register(bondn, bond_after)

end