function [ classseq ] = classseq( seq, lookup )
%classseq( seq, lookup ) create class sequence from aa sequence "seq" and fragment class
%lookup "lookup"

numit = size(lookup,5); % number of iterations
numaa = length(seq);

classseq = zeros(numaa-3,numit);

for i=1:numaa-3,
    aa1 = seq(i);
    aa2 = seq(i+1);
    aa3 = seq(i+2);
    aa4 = seq(i+3);

    if aa1 > 20 || aa2 > 20 || aa3 > 20 || aa4 > 20,
        continue;
    end
    
    classseq(i,:) = lookup(aa1,aa2,aa3,aa4,:);
    
end