function [ bondr, rotm ] = bond_normalizerotation( bond )
%[ bondr, rotm ] = bond_normalizerotation( bond ) 
%    bond == bondr*rotm'
%    bond*rotm == bondr
    
    ex = [1;0;0];
    r1 = bond(1,:)'/norm(bond(1,:));
    r2 = bond(2,:)'/norm(bond(2,:));

    % rotate so r1 || ex
    crx = cross(r1,ex); % = r1 x ex
    rotx = geom.rotmat(crx, acos(r1'*ex)); % rotate around cross product by angle given by dot product    
        
%     roty = geom.rotmat(ex, acos(r13rotyz'*r13roty/(norm(r13roty)*norm(r13rotyz))));
    
    r2x = rotx*r2;   % rotate r2 by rotx    
    
    % rotate rotated r2 around x axis so it lies in the xy plane    
    roty = geom.rotmat(ex, -atan(r2x(3)/r2x(2)));

    rotm = (roty*rotx)';

    bondr = bond*rotm;
end

