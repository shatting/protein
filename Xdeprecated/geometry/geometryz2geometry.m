function [ geom] = geometryz2geometry( geomz )
%GEOMETRYZ2GEOMETRY inverse of geometry2geometryz

nfrag = length(geomz.co)-1;

geom.c = int8(round(geomz.co(1:end-1)*100));
geom.co = geomz.co(1:end-1);

geom.cp = int8(round(geomz.co(2:end)*100));
geom.cpo = geomz.co(2:end);

y = sys10toX(geomz.s,[3 3 3]);
for i = 1:nfrag
      alphao(i,1) = geomz.alphaso(y(i,2) + (y(i,1) - 1)*3, y(i,3) + 2);
end

geom.alphao = alphao;

alpha = geomz.alphaso;
cz = [alpha(:,1:2),cos(alpha(:,3:5))];
sz = [alpha(:,1:2),sin(alpha(:,3:5))];
si = zeros(size(geomz.s));
ci = zeros(size(geomz.s));
sizes = size(geomz.s,2);
for i = 1:sizes
     si(i)  = sz(y(i,2) + (y(i,1) - 1)*3,y(i,3) + 2);
     ci(i) =  cz(y(i,2) + (y(i,1) - 1)*3,y(i,3) + 2);
end
% this is where I think we're supposed to be using the efficient formulas for to,
% tpo:
ti = zeros(sizes,1);
tprime = zeros(sizes,1);
for i = 1:sizes;
  tprime(i) = (ci(i)*(1 - geomz.zo(i)^2) - 2*si(i)*geomz.zo(i))/(1 + (geomz.zo(i))^2);
  ti(i) = (si(i)*(1 - geomz.zo(i)^2) + 2*ci(i)*geomz.zo(i))/(1 + (geomz.zo(i))^2);
end

geom.t = int8(round(100*ti));
geom.to = ti;
geom.tp = int8(round(100*tprime));
geom.tpo = tprime;

