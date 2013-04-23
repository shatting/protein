%check:
loadgeomdb;

problems = zeros(1,geomdb.nch); %rotation false
problems2 = zeros(1,geomdb.nch); % class sequence different
rmsd1 = zeros(1,geomdb.nch);
rmsd2 = rmsd1;
for i = 1:geomdb.nch  
    gd = geomdb.geomdata{i};
    s1 = HEL27_classifier(gd); 
    
    bond = double(gd.bond)./repmat(sqrt(sum((gd.bond).^2,2)),1,3);
    bondZ = bond_rotate(bond);
    if abs(bondZ(2,3)) > .00001
     problems(i) = 1;
    end
    
    [Q,q,rmsd1(i)] = register(bond,bondZ);
    rmsd2(i) = rmsd(bond,bondZ);
    
    s2 = HEL27_classifier(bond2geomstruct(bondZ));
    problems2(i) = sum(s1~=s2);    
end
probrot = find(problems ~= 0)
probclass = find(problems2 ~= 0)
meanrmsd1 = mean(rmsd1)
meanrmsd2 = mean(rmsd2)
