function s = mat2mol2(mat, title)

%header
s = ['MOLECULE ',title];
s = sprintf('%s\n%s\n',s,mat.header);

nat = length(mat.atomr);

%atoms
for at = 1:nat,
    a = mat.atomr(at);
    
    if isfield(a,'altloc'),
        s = [s sprintf('ATOM  %5i%5s%1s%3s %1s %4i%i   %8.3f%8.3f%8.3f%6.2f%6.2f%s\n',a.i,a.name,a.altloc,a.aaname,a.ch,a.aanum,a.icode,a.x,a.y,a.z,a.occ,a.bval,a.footn)];
    else % those from our data
        s = [s sprintf('ATOM  %5i%5s %3s  %4i    %8.3f%8.3f%8.3f%6.2f\n',a.i,a.name,a.aaname,a.aanum,a.x,a.y,a.z,a.occ)];
    end
end

s = [s sprintf('END\n')];

end