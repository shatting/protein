function [sane] = geom_test_suite(detailed,optidata)

if (nargin<1)
    detailed = 0;
end

sane = ~geom.bond2coords_test(detailed);
sane = sane && ~geom.bond_normalizerotation_test(detailed); 
sane = sane && ~geom.coords_register_test(detailed);
sane = sane && geom.ztransform_test(detailed);
sane = sane && protein.HEFG.HEFGSClassifier(1).test;
sane = sane && protein.HEFG.HEFGSClassifier(0).test;
sane = sane && geom.Coords.test_roundtrip(detailed);
if (nargin<2)
    dprintf('seq_pots_test skipped as no optidata passed to function');
else
    sane = sane && opti.seq_pots_test(optidata,detailed);
end

dprintf('- geom_test_suite.m');
if sane,    
    dprintf('PASSED.');
else
    dprintf('PROBLEMS FOUND!');
end