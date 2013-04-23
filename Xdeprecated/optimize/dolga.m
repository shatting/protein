function [GDT,LGArmse,lgacoords] = dolga(pre,datap)
% datap = data{p}
% pre = prediction in bonds

cals = calphas(pre);
mat1 = bond2mat(datap);
mat2 = mat1;
for i = 1:size(cals,1)
    mat2.atomr(i).x = cals(i,1);
    mat2.atomr(i).y = cals(i,2);
    mat2.atomr(i).z = cals(i,3);
end

% get into the correct directory
% I have no idea how to do this properly
% see if I am at home:
global proteinroot;
athome = 0;
if min(proteinroot == 'C:\Users\Owner\prot')
    athome = 1;
end
if athome == 1
    cd C:\Users\Owner\prot\ginny
end

filename = writelgainput(mat1,'template',mat2,'model','firsttry');

diary('lgaoutputb');
runlga(filename,'-3 -o2'); %options are supposed to come in the ''. want '-3' (and maybe -atom:CA -o1', but I think they're the defaults)
diary off
% get GDT
filename = 'lgaoutputb';
[GDT,LGArmse] = readoutGDT(filename);
% gdt LGA's prediction of best superposition of template onto model

cd ..
movefile('stephan\lga\lga_bin\firsttry.pdb','ginny\firsttry')
cd ginny
filename = 'firsttry';
lgacoords = readoutLGA(filename);
fclose('all');
delete 'lgaoutputb';
delete 'firsttry'
