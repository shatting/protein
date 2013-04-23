function classavgs = classaverages(newclasses,nums,data)
% finds average value of c, cp, z for each of the nums classes


ncl =  nums; %number of classes - 24
avgc = zeros(1,ncl); 
avgcplus = zeros(1,ncl);
avgz = zeros(1,ncl);
numg = zeros(1,ncl);

global alphas;

for p = 1:size(data,2)
   if mod(p,100) == 0
       p
   end
   if newclasses == 0
       geo = bond2geometry(data{p}.bond);
       geo = geometry2geometryz(geo, alphas);
       geo.cpo = geo.co(2:end);
       geo.co = geo.co(1:end-1);
       claseq = seq2sclasses(data{p}.bond);
       geo.z = geo.zo;
   else
       geo = bond2geo(data{p}.bond);
       claseq = geo.tcl; 
       geo.z = geo.tqu;
   end

   sizeseq = size(claseq,1); %should be 3 less than size sequence, bc 4-fragments
   
   for i = 1:sizeseq
       g = claseq(i);
       if ((g > 0) && (g < ncl + 1)) %just in case!
         avgc(g) = avgc(g)+geo.co(i);
         avgcplus(g) = avgcplus(g)+geo.cpo(i);
         avgz(g) = avgz(g)+geo.z(i);
         numg(g) = numg(g) + 1;
       end
   end %ends for i
   
end %ends for p 

%now just find averages
for g = 1:ncl
    if numg(g) > 0
    avgc(g) = avgc(g)/numg(g);
    avgcplus(g) = avgcplus(g)/numg(g);
    avgz(g) = avgz(g)/numg(g);
    end
end

classavgs = [avgc',avgcplus',avgz'];
Readmeavgs = ['columns describe avg c, cp, z geometry values for all classes'];

save classavg classavgs Readmeavgs


