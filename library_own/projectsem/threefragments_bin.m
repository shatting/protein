
load RT

tic
%sample = sampledata(data,10000);
sample = data;

nbins = 10;

angle = zeros(24,24,24, nbins);

numpr = length(sample);

for pr=1:numpr,
    
    if (mod(pr,100) == 0)
        disp(pr);
    end
    
    prcell = sample{pr};
    
    seq = prcell.seq;
    bond = double(prcell.bond);
    angles = getangles(bond);
    
    numaa = length(seq);    

    for i=1:numaa-2,
        aa1 = seq(i);
        aa2 = seq(i+1);
        aa3 = seq(i+2);
        
        mina = min(angles);
        maxa = max(angles);
        range = maxa-mina;
        
        bin = ceil(((angles(i)-mina)/range)*nbins);
        
        bin = min(max(bin,1),10);
        
        angle(aa1,aa2,aa3,bin) = angle(aa1,aa2,aa3,bin) + 1;
    end
    
end

count = countfrags(angle);

t = toc
toc
disp(sprintf('proteine/s: %f',numpr/t));
count