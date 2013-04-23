function [ num, midx ] = numsclasseslengthn( data, new_deal_data, sclasses )
%NUMSCLASSESLENGTHN 
% num(n,s) = maximum consecutive occurence length of HEL27 class sets sclasses{s} in proteins of
% length n
% i.e. sclasses = {[1,2,3],[4,5,6]} or sclasses = {1} or sclasses = {[1,2,3,4,5]}
% or numsclasseslengthn( data, new_deal_data, {[1,4,7,10,13,16,19,22,25]} )

nc = length(data);
nsets = length(sclasses);

dprintf('getting max sequence lengths from %i chains.',nc);
dprintf('number of class sets: %i',nsets);
for i=1:nsets, tightmat(sclasses{i}); end

lengths = zeros(1,nc);
for i=1:nc,
    lengths(i) = length(data{i}.seq);
end
maxl = max(lengths);


num = zeros(maxl,nsets);
midx = num;
figure(99);plot(0,0);hold on;
for i=1:nc,
    
    if (~mod(i,100)), dprintf('chain %i',i); end
    
    seq = data{i}.seq;
    l = length(seq);    
    
    %gamma = seq2gammaclasses(seq, new_deal_data);
    gamma = seq2sclasses(seq, data{i}.bond);
    for s=1:nsets,        
        n = maxsubseq(gamma,sclasses{s});
        if (num(l,s)<n)
            num(l,s) = n;
            midx(l,s) = i;
        end
    end
    
end


end
