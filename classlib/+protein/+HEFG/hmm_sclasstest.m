%options
nseqs = 100;
[chsidx, chs] = data.GeomDB.db.getsizencell(30,500);

dprintf('Using data analysis potential:');
new_deal_data_info_print(new_deal_data);

gclass = class.HEFGGammaClassifier(new_deal_data);
sclass = class.HEFGSClassifier(new_deal_data.diagbds);

HMM = hmm.HMM(sclass,gclass,data.FeatureDB.db);
distm = sclass.getdistm;
nchs = length(chs);

numhmmerr = 0;

resultdata = struct;
resultdata.hmmtime = 0;
resultdata.disttime = 0;
resultdata.naa = 0;
resultdata.dist = [];
resultdata.skipped = 0;
resultdata.chsidx = 0;
resultdata.avgdist = 0;
resultdata.avgdistnorm = 0;
resultdata.truesseq = 0;
resultdata.truesseqfound = 0;
resultdata.HMMSeqs = [];
resultdata.minnormdist = 0;
resultdata.mindist = 0;
resultdata(nchs).naa = 0;

skip = zeros(1,nchs);

parfor i=1:nchs,
%     if (~mod(i,round(nchs/10)))
%         dprintf('chain %i of %i (%i%%), skipped %i.',i,length(chs),round(i/nchs*100),sum(skip));
%     end     
    
    if (any(chs(i).seq > 20))  
        resultdata(i).skipped = 1;
        continue;
    end
        
    tic
    try
        HMMSseqs = HMM.getnbest(chs(i),nseqs);   
    catch e
        dprintf('something wrong with chain #%i.',chsidx(i));
        resultdata(i).skipped = 1;
        numhmmerr = numhmmerr +1;
        continue;
    end
    resultdata(i).skipped = 0;
    resultdata(i).hmmtime = toc;
    
    %dprintf('chain %i of %i took %s',i,length(chs),showtime(resultdata(i).hmmtime));
    resultdata(i).naa = length(chs(i).seq);
    resultdata(i).chsidx = chsidx(i);
    
    trueSseq = sclass.classify(chs(i));
    resultdata(i).truesseq = trueSseq;
    tic;
    dists = zeros(1,nseqs);
    for j=1:nseqs,
        dists(j) = hmm_sclassseqdist(HMMSseqs(j,:)',trueSseq,distm);            
    end
    if (any(dists==0))
        resultdata(i).truesseqfound = 1;
    else
        resultdata(i).truesseqfound = 0;
    end
    resultdata(i).dist = dists;
    resultdata(i).disttime = toc;
    resultdata(i).HMMSeqs = HMMSseqs;
    
    resultdata(i).avgdist = sum(dists)/nseqs;
    resultdata(i).avgdistnorm = resultdata(i).avgdist/(resultdata(i).naa-3);
    resultdata(i).minnormdist = min(dists)/(resultdata(i).naa-3);
    resultdata(i).mindist = min(dists);
%     if naa(i) < 20 & 0,
%         dprintf('chain #%i with %i naa. (first row: true s seq, other rows are predictions. first column: distance of prediction)',chsidx(i),naa(i));
%         disp([[0,-1,trueSseq'];zeros(1,naa(i)-1)-1;[dist(i,:)',zeros(nseqs,1)-1,HMMSseqs]]);            
%     end
end
% 
skipped = [resultdata.skipped];
resultdata = resultdata(~skipped);

dprintf('done.\nskipped %i.',sum(skipped));
dprintf('total time for HMM calculation: %s',showtime(sum([resultdata.hmmtime])));
dprintf('total time for dist calculation: %s',showtime(sum([resultdata.disttime])));
dprintf('HMM threw %i exceptions.',numhmmerr);
truesseqfound = [resultdata.truesseqfound];
dprintf('True sseq found in %i chains (%.2f%% of cases).',numhmmerr,sum(truesseqfound)/length(resultdata)*100);

avgdistnorm = [resultdata.avgdistnorm];
[m,mi] = min(avgdistnorm);
minres = resultdata(mi);
dprintf('minimum normalized average distance found: %.2f (chain #%i, %i naa)',m,minres.chsidx,minres.naa);
%dprintf('(first row: true s seq, other rows are predictions. first column: distance of prediction)');
%disp([[0,-1,minres.truesseq'];zeros(1,minres.naa-1)-1;[minres.dist',zeros(nseqs,1)-1,minres.HMMSeqs]]);   
minnormdist = [resultdata.minnormdist];
[m,mi] = min(minnormdist);
minres = resultdata(mi);
[~,minhmm] = min(minres.dist);
dprintf('minimum normalized distance found: %.2f (chain #%i, %i naa, HMM #%i)',m,minres.chsidx,minres.naa,minhmm);
dprintf('(first row: true s seq, other rows are predictions. first column: distance of prediction)');
disp([[0,-1,minres.truesseq'];zeros(1,minres.naa-1)-1;[minres.dist',zeros(nseqs,1)-1,minres.HMMSeqs]]);   

figure
naa = [resultdata.naa];
[a,b] = sort(naa);
subplot(3,1,1);
plot(naa(b),[resultdata(b).avgdist],'.k');
hold on;
plot(naa(b),[resultdata(b).mindist],'.r');
title('average (min) sseq distance versus chain length');
hold off;

subplot(3,1,2);
plot(naa(b),avgdistnorm(b),'.k');
hold on;
plot(naa(b),minnormdist(b),'.r');
title('average (min) sseq distance divided by chain length versus chain length');

subplot(3,1,3);
plot(naa(b),minnormdist(b)./avgdistnorm(b),'.g');
title('minimum/average sseq distance divided by chain length versus chain length');
