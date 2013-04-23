function [ pot ] = suffs2pots( suffdata, cl )
%GETPOT Summary of this function goes here
%   Detailed explanation goes here

suffgrps = suffdata.suffs;
ncl = max(cl);

for c=1:ncl,
    clind = find(cl==c);
    
    clsuff = suffgrps{clind(1)};
    if (length(clind)>1)
        for i=clind(2:end)',
            addsuff = suffgrps{i};
            clsuff=suffadd(clsuff,addsuff);
        end
    end
    classsuffs(c) = clsuff;
    pot(c) = suffpot(clsuff);
end


end
