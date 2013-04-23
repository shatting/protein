function [ probsfound, rms ] = coords_register_test(detailed)

dprintf('- coords_register.m: ');

for i=1:data.GeomDB.nch,
    coords = data.GeomDB.db.chains{i}.coords;
    shiftvec1 = 25*(rand(1,3)-0.5);
    shiftvec2 = 25*(rand(1,3)-0.5);
    
    bond = geom.coords2bond(coords); %    
    % get normalized bonds
    [bondn, rotm] = geom.bond_normalizerotation(bond);    
    
    sncoords = geom.bond2coords(bondn,shiftvec1); % shifted normalized coords
    socoords = geom.bond2coords(bond,shiftvec2);  % shifted original coords
        
    % find rotation and translation that normalizes shifted original coords
    % by registering
    [ Q, q, rrms(i) ] = geom.coords_register( sncoords, socoords );
    
    %rotdiff(i) = sum((rotm(:)- Q(:)).^2); % frobenius norm of rotation matrix difference     
    
    rms(i) = geom.rmsd(sncoords, geom.shift(socoords*Q',q));
    naa(i) = length(data.GeomDB.db.chains{i}.seq);
end

probsfound = max(rms) > 1e-8;

if (nargin > 0 && detailed) || probsfound
    [a,b] = sort(naa);

    figure;    
    plot(rms(b)/max(rms),'*b');
    hold on;
    plot(rrms(b)/max(rrms),'*k');
    plot(naa(b)/max(naa),'+r');
    title('coords_register.m');
    legend({sprintf('rms (max = %g)',max(rms)),sprintf('reported rms (max = %g)',max(rrms)),sprintf('naa (max=%i)',max(naa))});
    
    dprintf(' * maximum rmsd: %g',max(rms));
    dprintf(' * maximum difference rmsd-reported rmsd: %g',max(abs(rms-rrms)));    
end

if (probsfound)
    dprintf('PROBLEMS FOUND!');
else
    dprintf('PASSED.');
end