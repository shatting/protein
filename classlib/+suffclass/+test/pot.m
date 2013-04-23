% SUFF_POTTEST Test suff_data2cl and suff_data2v and compare classifiaction results.
n=100000;
ncat = 10;
nnum = 10;
ncl = 10;
numcenters = (rand(nnum,1)-0.5)*10.*rand(10,1);
cattype = randperm(ncat) + 10;

ftype = [ zeros(1,nnum) cattype ];

%num
data=rand(n,nnum) + repmat(numcenters',n,1);
%cat
datac=ceil(rand(n,nnum).*repmat(cattype,n,1));

data=[data datac];

cl = ceil(rand(n,1)*ncl);

suff=suffstat(ncl,ftype);

suff=suffstat(suff,data,n,cl);

pot = suff_pot(suff);

cpred = suff_data2cl( data, pot );
cpred2 = suff_data2v( data, pot, [pot.numpot.freq]/sum([pot.numpot.freq]) );

disp('target <-> suff_data2cl');
getfreq([cl cpred])

disp('target <-> suff_data2v');
getfreq([cl cpred2])

disp('suff_data2cl <-> suff_data2v');
getfreq([cpred cpred2])