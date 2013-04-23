function [ linear,seq ] = fragdb2linear( fragdb )
%FRAGDB2LINEAR Summary of this function goes here
%   Detailed explanation goes here

nfrag = fragdb.nfrag;
nfeat = fragdb.nfeat;
npri = prod(size(fragdb.frags));

linear = int16(zeros(nfrag,nfeat));
seq = uint8(zeros(nfrag,4));

c = 1;
for aa1=1:20,
    for aa2=1:20,
        for aa3=1:20,
            for aa4=1:20,                
                frags = fragdb.frags{aa1,aa2,aa3,aa4};
                if (isstruct(frags)),
                    t = c + frags.count;
                    linear(c:t-1,:) = int16(frags.data);
                    seq(c:t-1,:) = uint8(repmat([aa1,aa2,aa3,aa4],t-c,1));
                    c = t;
                end
            end
        end
    end
end