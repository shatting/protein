function [ mat ] = bond2mat( data )
%BOND2MAT Summary of this function goes here
%   Detailed explanation goes here

%coords = bond2coords(rotatematrix(double(data.bond)/10));
%coords = rotatematrix(double(data.bond)/10);
bond = rotatematrix(data.bond);
coords = calphas(bond);

naa = length(data.seq)% - 1;
mat = struct;
mat.header = sprintf('HEADER chain %s, only c-alphas',data.name);

for i=1:naa
    % atom record
    atomr = struct;     

    %s = [s sprintf('ATOM  %5i%5s%1s%3s %1s %4i%i   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n',a.i,a.name,a.altloc,a.aaname,a.ch,a.aanum,a.icode,a.x,a.y,a.z,a.occ,a.bval,a.footn)];    
    atomr.i = i;
    atomr.name = 'CA';
    atomr.aaname = aaname(data.seq(i));
    atomr.ch = 1;
    atomr.aanum = i;
    atomr.x = coords(i,1);
    atomr.y = coords(i,2);
    atomr.z = coords(i,3);    
    atomr.occ = 1;
    
    % save atom record
    mat.atomr(i) = atomr;
end


end