if exist('data')~=1, 
  tic
  clear
  load RT.mat 
  disp('data loaded')
  showtime
end;

tic


%sample = sampledata(data,10000);
sample = data;

nbin=10;
histang = zeros(20,20,20,nbin);

numpr = length(sample);
count=0;

for pr=1:numpr,
    
    if (mod(pr,100) == 0)
        disp(pr);
    end
    
    prcell = sample{pr};
    
    seq = prcell.seq;
    bond = double(prcell.bond);
    angles = getangles(bond);
    angmin=min(angles);
    angmax=max(angles)*1.00001;
    binlen=(angmax-angmin)/nbin;
    
    numaa = length(seq);    

    for i=1:numaa-2,
        aa1 = seq(i);
        aa2 = seq(i+1);
        aa3 = seq(i+2);
                
        bin=round((angles(i)-angmin)/binlen+0.5);
        if  histang(aa1,aa2,aa3,bin)==0, count=count+1; end;
        histang(aa1,aa2,aa3,bin) = histang(aa1,aa2,aa3,bin)+1 ;
    end
    
end

t = toc
showtime
disp(sprintf('proteine/s: %f',numpr/t));
count
