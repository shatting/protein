function [threefrags,Vfrag] = get3fragpot(p,data,pair2cl,potc)
% calculates the value of the 3-frag potential for a protein p
%tf has tf.pair2cl and tf.potc (tf = threefrags_interface...)
% bc suffpot defines pot.L as inv(chol(pot.cov)'), I think L is what I need
% to use

%pair2cl = tf.pair2cl;
wnormmax = 91125;
r2max = 7^2;

%find all 3-fragments and there classes
%bond = double(data{p}.bond)/10; %/10 might be important for ampl! 
bond = rotatematrix(data{p}.bond); %should have the same 3frag value as above, but just so everything matches
ca = bond2coords(bond);   
gramidx = [6,11,12,16,17,18,21,22,23,24];
aseq = double(data{p}.seq);
naa = size(aseq,1);
dmin = 1;
f = 0;
threefrags = zeros(naa^2,16); %first two columns: indices of frag, next two: aa types of frag, fifth: class of frag. next 11: feats
%rl = zeros(naa^2,1);

%for j = (3 + dmin + 1) : naa - 2 %minus 2, because only naa-3 bonds
%    for i = 1: (naa - 3 - dmin)
%        if (j - i) >= (3 + dmin)
ninter = 1;

for i = 1:(naa - 5 - ninter)
    for j = (i + 3 + ninter):naa - 2
            if ((aseq(i) < 21) && (aseq(j) < 21))
              f = f + 1;
 %             i,j
              threefrags(f,1:5) = [i,j,aseq(i),aseq(j),pair2cl(aseq(i),aseq(j))];
              
              %these formulas are taken from getrawdata.m !!
              
              r = ca(j+1,:) - ca(i+1,:); %% = bond(i+1,j+1)
              
              nr = norm(r);
%              rl(f) = nr;
              
              v = [bond(i:i+1,:); r; bond(j:j+1,:)];
              gram = v*v';
              feats = [gram(gramidx) nr]; %nr is gram(13)
              % check with ampl version of feats!!
%               feats(1), bond(i,:)*bond(i+1,:)'
%               feats(2), bond(i,:)*r'
%               feats(3), bond(i+1,:)*r'
%               feats(4), bond(i,:)*bond(j,:)'
%               feats(5), bond(i+1,:)*bond(j,:)'
%               feats(6), r*bond(j,:)'
%               feats(7), bond(i,:)*bond(j+1,:)'
%               feats(8), bond(i+1,:)*bond(j+1,:)'
%               feats(9), r*bond(j+1,:)'
%               feats(10), bond(j,:)*bond(j+1,:)'
%               feats(11), sqrt(r*r')
              
              
              
              threefrags(f,6:16) = feats;
            end
    %end
    end
end

f

num3frags = 0;
for i = 1:(naa - 5 - ninter)
    num3frags = num3frags + ((naa - 2) - (i + 3 + ninter) + 1);
end
num3frags

threefrags = threefrags(1:f,:);
%rl = rl(1:f);

%find potential!
pot = potc; %tf.potc;
Vfrag = 0;

    % run through all frags
    %r2max = 7^2; %not sure about this, but took from wtfn function, and rawdata_r_2_to_7;
        for i = 1:f; 
            L = pot(threefrags(i,5)).L;
            me = pot(threefrags(i,5)).mean;
            ent = pot(threefrags(i,5)).ent;
            if r2max > threefrags(i,16)^2  %(rl(i)^2)
                wl = (r2max -  threefrags(i,16)^2)^3;
            else
                wl = 0;
            end
            wlnorm = (wnormmax - wl)/wnormmax;
            Vfrag = Vfrag + wlnorm*(norm(L'*(threefrags(i,6:16)' - me),2)^2 + ent); %not sure if plus or minus ent
            %Vfrag = Vfrag + norm(L'*(threefrags(i,6:16)' - me),2) + ent;
        end % end i
end 
    
