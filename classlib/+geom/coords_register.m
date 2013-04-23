function [ Q, q, rmsd ] = coords_register( coordsx, coordsy )
%REGISTER find optimal superposition of coordsy on coordsx wrt. RMSD using the prokrustes
%method ie. kabsch' algorithm
%   finds optimal transform s.t. coordsx ~ coordsy*Q + q, or in other words
%                                rmsd(coordsx,shift(rotate(coordsy,Q),q))=min!
% INPUT:
%   coordsx(i,:)  ... ith original vector
%   coordsy(i,:)  ... ith to be superpositioned vector
% OUTPUT:
%   Q, q    ... x ~ shift(rotate(y,Q),q)
%   rmsd    ... root mean squared deviation

if(strcmp(class(coordsx),'geom.Coords') || strcmp(class(coordsx),'geom.Chain')),
    coordsx = coordsx.coords;
end
if(strcmp(class(coordsy),'geom.Coords') || strcmp(class(coordsy),'geom.Chain')),
    coordsy = coordsy.coords;
end

xmean = mean(coordsx);
ymean = mean(coordsy);

n = size(coordsx,1);

A = zeros(size(coordsx,2));
B = 0;

for i=1:n,
    A = A + (coordsx(i,:)-xmean)'*(coordsy(i,:)-ymean);
    B = B + sum((coordsx(i,:)-xmean).^2) + sum((coordsy(i,:)-ymean).^2);
end

[U, S, V] = svd(A);

%we dont want a reflection
%(http://cnx.org/content/m11608/latest/#MatrixAlignment)
d = sign(det(A)); 
D = eye(3); 
D(3,3) = d;

Q = U*D*V';
q = xmean - ymean*Q';
rmsd = norm(sqrt(1/n*(B - 2*trace(D*S)))); % norm(.) to get rid of possible tiny imaginary values