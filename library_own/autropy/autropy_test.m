clear autropy_test*

n=1000;
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

autropy_test_common