%SEQ2FRAGSEQ Converts aa sequence seq(:) to its corresponding 4-fragment
%sequence.
% INPUT: 
%       seq(i)          ith amino acid in the sequence
% OUTPUT:
%       fragseq(i,:)    ith 4-fragment
function fragseq = seq2fragseq(seq)            
    n=length(seq);

    fragseq = zeros(n-3,4);

    for i=1:n-3,
        fragseq(i,:) = seq(i:i+3);
    end
end

