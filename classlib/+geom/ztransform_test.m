function sane = ztransform_test(detailed)

dprintf('- ztransform_test.m');

ds = data.FeatureDB.db.getdatasetclone({'tp','t','beta'});
ds.ztransform(protein.HEFG.HEFGSClassifier(1),1);
tic;
b = ds.betatransform(protein.HEFG.HEFGSClassifier(1));
ttp = [cos(b), sin(b)];
t1 = toc;
maxabs1 = max(abs(ds.getdata('beta') - b));

ds = data.FeatureDB.db.getdatasetclone({'tp','t','beta'});
ds.ztransform(protein.HEFG.HEFGSClassifier(0),1);
tic;
b = ds.betatransform(protein.HEFG.HEFGSClassifier(0));
ttp = [cos(b), sin(b)];
t4 = toc;
maxabs4 = max(abs(ds.getdata('beta') - b));

ds = data.FeatureDB.db.getdatasetclone({'tp','t','beta'});
ds.ztransform(protein.HEFG.HEFGSClassifier(0),1);
tic;
[t,tp] = ds.ttptransform_rational(protein.HEFG.HEFGSClassifier(0));
t2 = toc;
maxabs2t = max(abs(ds.getdata('t') - t));
maxabs2tp = max(abs(ds.getdata('tp') - tp));

ds = data.FeatureDB.db.getdatasetclone({'tp','t','beta'});
ds.ztransform(protein.HEFG.HEFGSClassifier(1),1);
tic;
[t,tp] = ds.ttptransform_rational(protein.HEFG.HEFGSClassifier(1));
t3 = toc;
maxabs3t = max(abs(ds.getdata('t') - t));
maxabs3tp = max(abs(ds.getdata('tp') - tp));

if nargin > 0 && detailed
    dprintf('[%s for ttp] maximum z/beta transform absolute deviation (diag): %d',showtime(t1),maxabs1);
    dprintf('[%s for ttp] maximum z/beta transform absolute deviation (orth): %d',showtime(t4),maxabs4);
    dprintf('[%s] maximum z/(t,tp) transform absolute deviation (diag): (%d,%d)',showtime(t3),maxabs3t,maxabs3tp);
    dprintf('[%s] maximum z/(t,tp) transform absolute deviation (orth): (%d,%d)',showtime(t2),maxabs2t,maxabs2tp);
end

sane = max([maxabs1,maxabs4,maxabs2t,maxabs2tp,maxabs3t,maxabs3tp]) < 1e-14;

if sane,    
    dprintf('PASSED.');
else
    dprintf('PROBLEMS FOUND!');
end

