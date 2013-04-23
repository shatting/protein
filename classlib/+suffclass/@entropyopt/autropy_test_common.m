%needs suff, cl, data

pot2 = suff_pot(suff);
pr = getfreq(cl)/length(cl);

options.use_entropy = 0;
options.ent_weight = 0;
[ cpred, potmin, potvals] = suff_potvals( pot2, data, pr, options );
options.ent_weight = 0.5;
[ cpred2, potmin, potvals] = suff_potvals( pot2, data, pr, options );
options.ent_weight = 1;
[ cpred3, potmin, potvals] = suff_potvals( pot2, data, pr, options );

[ Vgd ] = vgammadelta( pot2, data, cl )

%[ Vgdpred ] = vgammadelta( pot2, data, cpred )

[ ent, delta ] = getentropies( Vgd, pr )

options.use_entropy = 1;
options.ent_weight = 1;
[ cpred4, potmin, potvals] = suff_potvals( pot2, data, ent, options );

disp('entweight 0');
tightmat(getfreq([cl cpred]));
tightmat(-sort(-getfreq(cpred)));

disp('entweight 0.5');
tightmat(getfreq([cl cpred2]));
tightmat(-sort(-getfreq(cpred2)));

disp('entweight 1');
tightmat(getfreq([cl cpred3]));
tightmat(-sort(-getfreq(cpred3)));

disp('autropy');
tightmat(getfreq([cl cpred4]));
tightmat(-sort(-getfreq(cpred4)));