function [ D ] = kldivpairs( r )
%KLDIVPAIRS computes all pair divergences D(r(:,i)||r(:,j)) for i~=j
% D(i,j) = D(r(:,i)||r(:,j))

[nx,nr] = size(r);
D = zeros(nr);

for i=1:nr,
    q = r(:,i);
    idxp=setdiff(1:nr,i);
    pq = kldiv(r(:,idxp),q); % pq(1,i) = D(r(idxp(i))||q)
    D(idxp,i) = pq';
end


end
