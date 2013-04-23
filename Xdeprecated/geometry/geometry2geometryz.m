function [ geomz ] = geometry2geometryz( geom, alphas )
%GEOMETRY2GEOMETRYZ convert geometry with c, t, t' and phi to geometry with
%c, z
% vis http://de.wikipedia.org/wiki/Tangens_und_Kotangens#Rationale_Parametrisierung
% input: 
%    geometry struct with (at least) the fields
%        c: [453x1 int8]
%       co: [453x1 double]
%       cp: [453x1 int8]
%      cpo: [453x1 double]
%      phi: [453x1 int8]
%     phio: [453x1 double] ** used to find z!
%
%    alphas (see HEL27_definition.m, [-100,100])
% output:
%   struct with fields
%        co(1:nfrag+1): [-1,1] = cos(bondangle_i)  all cs! no distinction between c and cp
%        zo(1:nfrag)  : [-inf,inf] = tan((phi_i - alpha_i)/2)
%              alphaso: [alphas in radians [-pi,pi]]
%                    s: the s classes that were used. 
%                       we save s and alphas, for if
%                       the HEL classification changes, we will still be
%                       able to reconstruct the fold

    nfrag = length(geom.c);

    alphaso = [alphas(:,1:2),alphas(:,3:5)*pi/100];
    zo = zeros(nfrag,1);
    alphao = zeros(nfrag,1);

    s = HEL27_anal08([geom.c geom.cp geom.phi])';

    y = sys10toX(s,[3 3 3]);

    for i = 1:size(geom.c,1)
      alphao(i) = alphaso(y(i,2) + (y(i,1) - 1)*3, y(i,3) + 2);
      %zo(i) = tan(((geom.phio(i)-alphao(i))/2));
    end
    for i = 1:size(geom.c,1)
        zo(i) = tan((geom.phio(i) - alphao(i))/2);
    end
    %geomz.c = [geom.c;geom.cp(end)]; dont want to store that here since we
    %dont get them out of optimization
    geomz.co = [geom.co;geom.cpo(end)];
    %geomz.z; % there is no z since no fixed interval: tan -> (-inf,inf)
    geomz.zo = zo; 
    geomz.alphaso = alphaso;
    geomz.s = s;
    
end
