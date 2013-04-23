
% very simple example data. the additional entries in class 1 and 2 are for
% regularization and classification testing.
testdata = [repmat([1 2],5,1);[2 1]; [1 2]; repmat([2 1],5,1);[2 1]; [1 2];[1 1]; [2 2]];
testtarget = [repmat(1,5+2,1);repmat(2,5+2,1);repmat(3,2,1)];

suff = suffstat(3,[2 2]);
suff = suffstat(suff,testdata,size(testdata,1),testtarget);

suff.count

[ pclfeats,pcl,pfeats,pfeatscl,direct ] = suff_catprob(suff,[2 2])

pot = suff_pot(suff);

options.handleequalpotentials = 1;
[ cl, vmin, V, prob, Vdiff ] = suff_data2v([1 1; 1 2; 2 1; 2 2],pot,pot.catpot.freq/sum(pot.catpot.freq),options)

for f=1:size(testdata,1)
   pclfeats(:,f) = suff_catprob(suff,testdata(f,:));
end
