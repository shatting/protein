function [ geometry ] = bond2geometry( bond )
%[ geometry ] = bond2geometry( bond )
% geometry.c
% geometry.s
% geometry.cp
% geometry.sp
% geometry.t
% geometry.tprime

bond = double(bond);

numaa = size(bond,1) + 1;

bondn = bond./repmat(sqrt(sum(bond.^2,2)),1,3); %normalized bond vectors

Gn = createG(bondn);

geometry.c = Gn(1:end-2,2);
geometry.s = sqrt(1-geometry.c.^2);

geometry.cp = Gn(2:end-1,2);
geometry.sp = sqrt(1-geometry.cp.^2);

cplus=geometry.c.*geometry.cp; % c_i * c_{i+1}
%cplus = cplus(1:end-1);

splus=geometry.s.*geometry.sp; % s_i * s_{i+1}
%splus = splus(1:end-1);

geometry.t = (cplus-Gn(1:end-2,3))./splus;

%|a,b,c| = a*(bxc)
dets=zeros(numaa-3,3);
for i=1:numaa-3,
    a = bondn(i+1,:);
    b = bondn(i+2,:);    
    dets(i,1)=a(2)*b(3)-a(3)*b(2);
    dets(i,2)=a(3)*b(1)-a(1)*b(3);
    dets(i,3)=a(1)*b(2)-a(2)*b(1);
end;
dets = bondn(1:end-2,:)*dets';

geometry.tprime=diag(dets)./splus;


%test if correct
%c.^2+s.^2
%geometry.t.^2+geometry.tprime.^2