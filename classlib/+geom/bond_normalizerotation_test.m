function [ probsfound, rms, problems ] = bond_normalizerotation_test(detailed)
%BON_NORMALIZEROTATION_TEST Summary of this function goes here
%   Detailed explanation goes here

dprintf('- bond_normalizerotation.m: ');

problems = [];
for i=1:data.GeomDB.nch,
    % get a bond
    bond = data.GeomDB.db.chains{i}.bond;
    % normalize it
    [bondr, rotm] = geom.bond_normalizerotation(bond);
    if (norm(bondr(1,2:3)) > 1e-8 || norm(bondr(2,3)) > 1e-8 )
        problems(end+1) = i;
    end
    % check if rotation was correct
    rms(i) = geom.rmsd(geom.bond2coords(bond),geom.bond2coords(bondr*rotm'));
    naa(i) = length(data.GeomDB.db.chains{i}.seq);
end

probsfound = ~isempty(problems) || max(rms) > 1e-8;

if (nargin > 0 && detailed) || probsfound
    [a,b] = sort(naa);

    figure;
    plot(rms(b)/max(rms),'*b');
    hold on;
    plot(naa(b)/max(naa),'+r');
    title('bond_normalizerotation.m');
    legend({sprintf('rms (max = %g)',max(rms)),sprintf('naa (max=%i)',max(naa))});
    
    dprintf(' * maximum rmsd: %g',max(rms));
    dprintf(' * not in correct configuration: %i',length(problems));
end


if (probsfound)
    dprintf('PROBLEMS FOUND!');
else
    dprintf('PASSED.');
end

end

