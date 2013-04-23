% try to assess impact of setting bondlengths == 3.8a° on lga score

pdbfile = 'examplefiles\pdb110l.ent';

% convert pdb file to mat struct
origm = pdb2mat(pdbfile);
% set lengths to 3.8
modm = setbonds(origm,3.8);

% get mol2 strings for lga
origs = mat2mol2(origm,'110l.orig');
mods = mat2mol2(modm,'110l.mod');

disp('sum of angle differences before and after setbonds(.,3.8)');
bnorig = bond2geometry(mat2bond(origm));
bnmod = bond2geometry(mat2bond(modm));

cerror = sum(bnorig.c-bnmod.c)
terror = sum(bnorig.t-bnmod.t)
tprimeerror = sum(bnorig.tprime-bnmod.tprime)

modfile = '110lmod.110lorig';
fid = fopen(modfile,'w');
fprintf(fid,mods);
fprintf(fid,'\n');
fprintf(fid,origs);
fclose(fid);

movefile(modfile,'lga_bin/MOL2');

dprintf('run "./lga_bin/lga [params] %s > outfile"',modfile);