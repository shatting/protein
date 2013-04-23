function [probsfound, rms ]= bond2coords_test(detailed)

dprintf('- bond2coords/coords2bond:')

for i=1:data.GeomDB.nch,
    bond = data.GeomDB.db.chains{i}.bond;
    bond2 = geom.coords2bond(geom.bond2coords(bond,[-3 0 3]));
    rms(i) = geom.rmsd(bond,bond2);
end

probsfound = max(rms) > 1e-8;

if (nargin>1 && detailed)
    figure
    plot(rms);
    title('bond2coords/coords2bond');
    dprintf(' * maximum rmsd: %g',max(rms))
end

if (probsfound)    
    dprintf('PROBLEMS FOUND!');
else
    dprintf('PASSED.');
end