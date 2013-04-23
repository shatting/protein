
ch = data.GeomDB.db.chains{1};

d = ch.geometry.getdata({'c','cp','t','tp'});

c = [d(:,1);d(end,2)];
t = d(:,3);
tp = d(:,4);

bond = geom.geometry2bond(c,t,tp,ch.geometry.len,1);
bond2 = geom.geometry2bond(c,t,tp,ch.geometry.len);

[~,~,rmsd] = geom.coords_register(geom.bond2coords(bond.x),ch.coords)
%[x,y,rmsd2] = geom.coords_register(geom.bond2coords(bond2),ch.coords)