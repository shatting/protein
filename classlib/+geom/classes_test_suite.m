
nchains = data.GeomDB.nch;

dprintf('geom.Geometry and geom.Coords test suite.');
for i=1:nchains,
    if (~mod(i,round(nchains/10)))
       dprintf('%i%%',round(i/nchains*100));
    end
    
    problems(i,1) = geom.Coords.test_normalize(i);
    problems(i,2) = geom.Coords.test_registerto(i);
    [problems(i,3), x, Acond(i)] = geom.Coords.test_geomreconstruction(i);
    
    naa(i) = length(data.GeomDB.db.chains{i}.seq);
end

[a,b] = sort(naa);

nproblems = sum(problems(:));
if (nproblems > 0)
    dprintf('TEST FAILED! : %i problems',sum(problems(:)));
    dprintf('problems in normalization: %i',sum(problems(:,1)));
    dprintf('problems in registration: %i',sum(problems(:,2)));
    dprintf('problems in reconstruction: %i',sum(problems(:,3)));
else
    dprintf('test passed.');
end
