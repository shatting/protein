function [mi,ma] = hydrobylength(data)
% looking for some way of defining the boundaries for Vhy in opt. problem
% don't know if anything will come out here, but might as well try
% hydro constraint has more to do with aa composition, so maybe I should
% try something more complicated after this
% mi, ma are the coefficients of linear functions plotting relationship bw
% length and Vhy
% used in createproteinfile to define hydroup, hydrodown ampl params

load hydropots.mat

[lo,sh,ne]=sizefrags(data); %ne is number of frags longer than 0

avghpot = zeros(ne,1); %stores average Vhy value for each protein length
minhpot = inf(ne,1); %stores minimal Vhy value for each protein length
maxhpot = -inf(ne,1); %stores maximal Vhy value for each protein length
Stop = zeros(ne,1);

clf;
si = 0; %keeps track of total number of different fragment lengths

for n = sh:lo
     if rem(n,100)==0 %just to entertain!
        n
    end
    if numfragslengthn(n,data)== 0
        continue;
    end
         si=si+1; %si is the total number of fragment lengths (should be equal to ne, but just in case)
         Stop(si)=n; % keep track of the different lengths seen, for plotting
         [sizec,C]=fragmentslengthn(n,data); %finds vector C containing all fragments of this length
         for i = 1:sizec
             p = C(i);
             Vhydro = hydropotential(data,p,pot);
             avghpot(si) = avghpot(si) + Vhydro;
             if Vhydro < minhpot(si)
                 minhpot(si) = Vhydro;
             end
             if Vhydro > maxhpot(si)
                 maxhpot(si) = Vhydro;
             end
         end
         avghpot(si) = avghpot(si)/sizec;
end %for n = sh:lo

% get rid of extra zeros
avghpot = avghpot(1:si);
minhpot = minhpot(1:si);
maxhpot = maxhpot(1:si);
Stop = Stop(1:si);

plot(Stop,avghpot,'r');
hold on;
plot(Stop,minhpot,'g');
hold on;
plot(Stop,maxhpot,'b');

% find linear approximations
mi = polyfit(Stop,minhpot,1);
ma = polyfit(Stop,maxhpot,1);
minh = @(x) mi(1).*x + mi(2);
maxh = @(x) ma(1).*x + ma(2);

hold on;
fplot(minh,[sh lo],'b');
hold on;
fplot(maxh,[sh lo],'g');

         