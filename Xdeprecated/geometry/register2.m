function [ Q, q ] = register2( x, y )
%REGISTER2 minimiert = rmsd(x,shift(rotate(y,Q),q))

xmean = mean(x);
ymean = mean(y);

x = shift(x,xmean);
y = shift(y,ymean);

C = x'*y;

[V, S, W] = svd(C);

d = sign(det(C));
D = eye(3); D(3,3) = d;

Q = V*D*W';
q = (xmean - ymean*Q');