% compare angles
disp('angles before and after setbonds(.,3.8)');
bnorig = bond2geometry(mat2bond(setbonds(origm,3.8)));
bnb = bond2geometry(mat2bond(origm));

cerror = sum(bnorig.c-bnb.c)
terror = sum(bnorig.t-bnb.t)
tprimeerror = sum(bnorig.tprime-bnb.tprime)

% compare 
disp('angles after setbonds(.,3.8) and after 3.8*geometry2bondn()');
bnorig = bond2geometry(mat2bond(setbonds(origm,3.8)));
bnb = bond2geometry(3.8*geometry2bondn(bond2geometry(mat2bond(origm))));

cerror = sum(bnorig.c-bnb.c)
terror = sum(bnorig.t-bnb.t)
tprimeerror = sum(bnorig.tprime-bnb.tprime)