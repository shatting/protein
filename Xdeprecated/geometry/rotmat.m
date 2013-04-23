function [ m ] = rotmat( v, alpha )
%[ m ] = rotmat( v, alpha )
% m is the matrix that positively rotates around v by alpha

if (norm(v) < 1e-15)
	m = eye(3);
    return
end
    

v = v/norm(v);
c = cos(alpha);
cm = 1-c;
s = sin(alpha);
x = v(1);
y = v(2);
z = v(3);

% http://de.wikipedia.org/wiki/Rotationsmatrix
m(1,:) = [c + x^2*cm, x*y*cm - z*s, x*z*cm + y*s];
m(2,:) = [y*x*cm + z*s, c + y^2*cm, y*z*cm - x*s];
m(3,:) = [z*x*cm - y*s, z*y*cm + x*s, c + z^2*cm];