function [ n, i ] = maxsubseq( seq, subseqel )
% [ n, i ] = maxsubseq( seq, subseqel )
% n is length of longest subsequence of elements of subseqel within seq,
% and i the starting index of the last occurence of such a longest sequence
% ("last" because there could be more than 1..)

naa=length(seq);
n = 0;
i=1;
while i<=naa
    j=i;
    while (j<=naa && ismember(seq(j),subseqel)),
        j = j+1;        
    end
    if j-i>25,
      plot(naa,j-i+rand-1/2,'.');
    end;
    n = max(n, j-i);
    i=j+1;
end


end
