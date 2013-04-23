function [ mypoly, p ] = constr_hydro_getfitp(detailed)
%CONSTR_HYDRO_FINDBDS Summary of this function goes here
%   Detailed explanation goes here

chs = data.GeomDB.db.chains;
nch = length(chs);
sumbysize = zeros(1500,1);
numbysize = sumbysize;
for i=1:nch,
    coords = chs{i}.coords;    
    n = size(coords,1);
    if (n>1500)
        error('chain longer as assumed, please adjust constr_hydro_getfitp.m');
    end
    center = mean(coords,1);
    cvecs = coords - center(ones(n,1),:);
    distvar(i) = sum(sum(cvecs.^2,2))/n;
    sumbysize(n) = sumbysize(n) + distvar(i);
    numbysize(n) = numbysize(n) + 1;
    
end
meanbysize = sumbysize./numbysize;

datax = find(numbysize > 0);
datay = meanbysize(datax);

% now find interpolating polynom of order 2/3
p = polyfit(datax.^(1/3),datay,2);
mypoly = @(x) p(1)*x.^(2/3) + p(2)*x.^(1/3) + p(3);

if (nargin > 0 && detailed)
    plot(datax,datay,'.b');
    hold on;

    plot(datax,mypoly(datax),'r')
    hold off;
end

end

