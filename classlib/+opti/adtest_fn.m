function [ c, coords ] = adtest_fn( x, len, sseq, betabar, hydropotseq, chin, useintlabAD )
%ADTEST_FN 

nfrags = (length(x)-1)/2;

% get t, tp from z
[t,tp] = geom.ttptransform_rational(x(nfrags+2:end),sseq,cos(betabar'),sin(betabar'));

% get coords from c, t, tp
bond = geom.geometry2bond(x(1:nfrags+1),t,tp,len, useintlabAD);

coords = geom.bond2coords(bond);

% the constraints
cc = opti.constr_cc( coords );
hydro = opti.constr_hydro(coords,hydropotseq,chin);

c = [cc; hydro];

end

