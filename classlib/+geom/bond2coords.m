function [ coords ] = bond2coords(bond,x0)
%GETCOORDS Summary of this function goes here
%   Detailed explanation goes here

if (~isa(bond,'gradient')),
    bond = double(bond);
end

if(nargin < 2)
    x0 = [0 0 0];
end        

dograd = 0;
if (isa(bond,'gradient')),
    %testcoords = gradient(zeros(size(bond,1)+1,3));
    %testcoords(1,:) = gradient(x0);
    %testbond = bond;
    
    dograd = 1;
    nxgrad = size(bond.dx,3);
    naa = size(bond.x,1)+1;
    Dbond = bond.dx;    
    Dcoords = zeros(naa,3,nxgrad);
    bond = bond.x;
else
    naa = size(bond,1)+1;
end    

coords = zeros(size(bond,1)+1,3);
coords(1,:) = x0;

for i=2:naa,
    coords(i,:) = coords(i-1,:) + bond(i-1,:);    
    if (dograd)
        %testcoords(i,:) = testcoords(i-1,:) + testbond(i-1,:);
        Dcoords(i,:,:) = Dcoords(i-1,:,:) + Dbond(i-1,:,:);
    end
end

if (dograd)
    r.x = coords;
    r.dx = reshape(Dcoords,3*naa,nxgrad);
    coords = gradient(r);
end

end

