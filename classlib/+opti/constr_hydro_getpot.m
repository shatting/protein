function [ pot, chi ] = constr_hydro_getpot(detailed)
%CONSTR_HYDRO_GETPOT Summary of this function goes here
%   Detailed explanation goes here

[ chi ] = opti.constr_hydro_getfitp(detailed);

chs = data.GeomDB.db.chains;
nch = length(chs);

aadata = suffclass.suffstat(20,0); % only one feature (13.10.09)


for i=1:nch,    
    coords = chs{i}.coords;    
    n = size(coords,1);
    seq = chs{i}.seq;
    center = mean(coords,1);
    cvecs = coords - center(ones(n,1),:);
    idxok = seq<=20; % really would not need this since suffstat ignores all classes that have labels higher than max label at creation
    if (~sum(idxok))
        continue;
    end
    distaasquared = sum(cvecs(idxok).^2,2); % distance of aas to center squared
    
    aadata.adddata(distaasquared/chi(n),seq(idxok));
    
end
opts = struct;
opts.num_reg_weight = 0;
%opts.debug_usematcov = 1;
pot = suffclass.potential(aadata,opts);

end

