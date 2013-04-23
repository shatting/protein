function [bondZ] = rotatematrix(bond)
%rotates the bond matrix onto [[norm(bond1),0,0];[something something,0];...]
% important to make sure that the s sequence is the same after rotation!!!!
% see checkrotatematrix - it works!!!! and easy to explain mathematically!

bond = double(bond);
%normalizes the matrix first!! (not necessary, but convenient)
bond = (bond)./repmat(sqrt(sum((bond).^2,2)),1,3);

[th,ph,r]=cart2sph(bond(1,1),bond(1,2),bond(1,3));
t = -th;
Z = [[cos(t),sin(t),0];[-sin(t),cos(t),0];[0,0,1]];%matrix to rotate around Z axis
bondX = bond*Z;

[th,ph,r]=cart2sph(bondX(1,1),bondX(1,2),bondX(1,3));
t = ph;
Y = [[cos(t),0,-sin(t)];[0,1,0];[sin(t),0,cos(t)]]; % matrix to rotate around Y axis
bondY = bondX*Y;

t = atan(bondY(2,3)/bondY(2,2)); %chose this t, so that -bond(2,3)*sint + bond(3,3)*cost = 0 (just solve equation for t, and use tan = sin/cos)!
X = [[1 0 0];[0 cos(t) -sin(t)];[0 sin(t) cos(t)]]; %matrix to rotate around X axis 
bondZ = bondY*X;


