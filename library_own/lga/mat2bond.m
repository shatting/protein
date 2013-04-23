function [ bond, naa, aaidx, caidx ] = mat2bond( m )
%[ bond, naa, aaidx, caidx ] = mat2bond( m )
% 

aanum = [m.atomr(:).aanum];
naa = max(aanum);
caidx = zeros(0,1);
bond = zeros(naa-1,3);

for i=1:naa,
    idxaa = find(aanum==i);
    aaidx{i} = idxaa;
    natoms = length(idxaa);
    for a=1:natoms,
        if (strcmp(strtrim(m.atomr(idxaa(a)).name),'CA')),
            caidx(end+1) = idxaa(a);
        end
    end  
    if i>1,
        a1 = m.atomr(caidx(i-1));
        a2 = m.atomr(caidx(i));
        bond(i-1,:) = [a2.x-a1.x a2.y-a1.y a2.z-a1.z];
    end
end


end
