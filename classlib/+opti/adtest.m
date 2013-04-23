data.GeomDB.reloadch(7);
ch = data.GeomDB.db.chains{7};

cl = class.HEFGSClassifier(1);
ch.ztransform(cl,1);
sseq = ch.getdata({'clsz'});
betabar = cl.gettorsioncenters;
x0 = ch.getknitrofeatures;
xgrad = gradientinit(x0);

[pot,chi] = opti.constr_hydro_getpot(0);
hydropotseq = pot(ch.seq);
chin = chi(length(ch.seq));

fun = @(x,il)opti.adtest_fn(x,ch.geometry.len, sseq, betabar,hydropotseq, chin, il);

tic
dprintf('performing manual AD on function.');
[ cgrad, coordsgrad ] = fun( xgrad, 0 );
tim1 = toc;

tic
dprintf('performing intlab+manual AD on function.');
[ cgrad2, coordsgrad2 ] = fun( xgrad, 1 );
tim3 = toc;

tic
[ c, coords ] = fun( x0, 1 ) ;
tim2 = toc;

% test if correctly transformed + reconstructed
[~,~,recrmsd] = geom.coords_register(coordsgrad.x,coords);
[~,~,recrmsd2] = geom.coords_register(coordsgrad2.x,coords);

dprintf('performing numerical differentiation on function.');
[gradnum,err] = jacobianest(@(x)fun(x,0), x0);

dprintf('\nfunction evaluation took %.3fs with manual AD, %.5fs without (factor of %.2f)',tim1,tim2,tim1/tim2);
dprintf('manual AD reconstruction rmsd: %g',recrmsd);
dprintf('relative error in function value (cc,hydro)=(%g,%g)',abs(c(1)-cgrad(1).x)/c(1),abs(c(1)-cgrad(1).x)/c(1));
dprintf('relative error in gradient (cc,hydro)=(%g,%g)\n',norm(cgrad(1).dx-gradnum(1,:))/norm(gradnum),norm(cgrad(2).dx-gradnum(2,:))/norm(gradnum));

dprintf('function evaluation took %.3fs with manual+intlab AD, %.5fs without (factor of %.2f)',tim3,tim2,tim3/tim2);
dprintf('intlab AD reconstruction rmsd: %g',recrmsd2);
dprintf('relative error in function value (cc,hydro)=(%g,%g)',abs(c(1)-cgrad2(1).x)/c(1),abs(c(1)-cgrad2(1).x)/c(1));
dprintf('relative error in gradient (cc,hydro)=(%g,%g)',norm(cgrad2(1).dx-gradnum(1,:))/norm(gradnum),norm(cgrad2(2).dx-gradnum(2,:))/norm(gradnum));