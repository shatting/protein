function [ bond,rotm ] = bond_rotate( bond )
%[ bond,rotm ] = bond_rotate( bond ) rotate bonds to canonical
%configuration.

bond = double(bond);

ex = [1;0;0];
r1 = bond(1,:)'/norm(bond(1,:));
r2 = bond(2,:)'/norm(bond(1,:));

% rotate so r1 || ex
crx = crossmat(r1)*ex; % = r1 x ex
rotx = rotmat(crx, acos(r1'*ex)); % rotate around cross product by angle given by dot product

% rotate rotated r2 around x axis so it lies in the xy plane
r2yz = rotx*r2; r2yz(1) = 0;
r2y = rotx*r2; r2y([1,3]) = 0;
roty = rotmat(ex, acos(r2yz'*r2y/(norm(r2yz)*norm(r2y))));

rotm = (roty*rotx)';

for i=1:size(bond,1),
    bond(i,:) = bond(i,:)*rotm;
end