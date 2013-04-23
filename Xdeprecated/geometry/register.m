function [ Q, q, rmsd ] = register( x, y )
%REGISTER find optimal superposition of y on x wrt. RMSD using the prokrustes
%method ie. kabsch' algorithm
%   finds optimal transform s.t. x ~ yQ + q, or in other words
%                                rmsd(x,shift(rotate(y,Q),q))=min!
% INPUT:
%   x(i,:)  ... ith original vector
%   y(i,:)  ... ith to be superpositioned vector
% OUTPUT:
%   Q, q    ... x ~ shift(rotate(y,Q),q)
%   rmsd    ... root mean squared deviation

xmean = mean(x);
ymean = mean(y);

n = size(x,1);

A = zeros(size(x,2));
B = 0;

for i=1:n,
    A = A + (x(i,:)-xmean)'*(y(i,:)-ymean);
    B = B + sum((x(i,:)-xmean).^2) + sum((y(i,:)-ymean).^2);
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