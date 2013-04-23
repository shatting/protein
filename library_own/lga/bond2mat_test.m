
% 1crn is only in DSbond.mat
i = findchain(data,'1crn');
% get data
b = data{i};

% convert to mat struct
bmat = bond2mat(b);

% get mat struct from pdb file, only c-alphas
dsmat = pdb2mat('examplefiles\1crn.pdb',true);

filename = writelgainput(bmat,'from data{}',dsmat,'from pdb file','1crndata.1crnpdb');

runlga(filename,'');

% funny how bond2coords(b.bond/10) doesnt work here, but
% bond2coords(b.bond)/10 does..
% resolved..has to be double(b.bond)/10
[Q,q,rmsd] = register(bond2coords(mat2bond(dsmat)),bond2coords(double(b.bond)/10))