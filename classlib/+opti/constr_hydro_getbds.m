function [  ] = constr_hydro_getbds(  )
%CONSTR_HYDRO_GETBDS Summary of this function goes here
%   Detailed explanation goes here

[ pot, fitp ] = opti.constr_hydro_getpot(1);

chs = data.GeomDB.db.chains;
nch = length(chs);
constrvals = [];
len = [];
for i=1:nch,    
    coords = chs{i}.coords;    
    n = size(coords,1);
    seq = chs{i}.seq;
    idxok = seq<=20; 
    if (any(~idxok))
        continue;
    end    
    
    constrvals = [constrvals; opti.constr_hydro( coords, pot(seq), fitp(n) )];
    len = [len,length(chs{i}.seq)];
end
[~,b] = sort(len);
x = len(b);
y = constrvals(b);

figure;
plot(x,y,'.');

[c,d] = sort(y);
idxmin = d(round(0.05*length(y)));
idxmax = d(round(0.95*length(y)));

end

