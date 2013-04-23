function [ vhydro ] = constr_hydro( coords, hydropotseq, hydro_chin )
%CONSTR_HYDRO Summary of this function goes here
%   Detailed explanation goes here
n = size(coords,1);
center = mean(coords,1);

cvecs = coords - center(ones(n,1),:);

cminuscbardivchi = sum(cvecs.^2,2)/hydro_chin;

% works as .L and .mean are 1x1
Ls = [hydropotseq.L]';
mus = [hydropotseq.mean]';
ents = [hydropotseq.ent]';
vali = (Ls .* (cminuscbardivchi-mus)).^2 - ents;


vhydro = sum(vali);

end

