function [ ccmax, len, cc, polycoeff ] = constr_cc_findmax( coords )
%CONSTR_CC_FINDMAX Summary of this function goes here
%   Detailed explanation goes here
ccmax = 0;
time = 0;
n = 1;
if (nargin > 0)    
    if (isa(coords,'geom.Coords'))
        coords = coords.coords;
    end        
    tic
    ccmax = opti.constr_cc(coords);
    time = toc;
    naa = size(coords,1);
else	
    chs = data.GeomDB.db.chains;
    for i=1:length(chs),
        if (~mod(i,round(length(chs)/20))),
            dprintf('calculating cc values in database %i%%',round(i/length(chs)*100));
        end
        len(i) = length(chs{i}.seq);
        tic
        cc(i) = opti.constr_cc(chs{i}.coords);
        time = time+toc;        
    end
    n = length(chs);
    naa = data.GeomDB.db.naa;
    
    [a,b] = sort(len);
    len = len(b);
    cc = cc(b);

    plot(len,cc,'.b');
    hold on;
    [polycoeff,s] = polyfit(len,cc,3);
    plot(len,polyval(polycoeff,len),'r');
    hold off;
    figure
    plot(len,cc./polyval(polycoeff,len));
end

dprintf('%s, processed %i chains, %.2f aa/sec.',showtime(time),n,naa/time);

end

