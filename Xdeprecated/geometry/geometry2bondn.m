function [ bond ] = geometry2bondn( geometry )
%GEOMETRY2BOND calculate normalized bond vectors from geometry c,t,t'
%   first bond vector is parallel to x-axis, second lies in xy-plane

if isfield(geometry,'z'),
    error('please use geometryz2geometry.m first.');
end

c = [geometry.co;geometry.cpo(end)];

s = sqrt(1-c.^2);

bond = zeros(size(c,1),3);

bond(1,1) = 1;
bond(2,1) = c(1);
bond(2,2) = s(1);
for i=1:length(c)-1,
    a(1,:) = bond(i,:)*crossmat(bond(i+1,:)'); % r_i x r_{i+1}
% var a11(i in ILIST):= bond2(i)*bond3(i+1)
% var a12(i in ILIST):=  
    a(2,:) = bond(i,:);
    a(3,:) = bond(i+1,:);
    b = [geometry.tpo(i)*s(i)*s(i+1); c(i)*c(i+1)+geometry.to(i)*s(i)*s(i+1); c(i+1)];
    bond(i+2,:) = (a\b)';
    bond(i+2,:) = bond(i+2,:)/norm(bond(i+2,:));
end
