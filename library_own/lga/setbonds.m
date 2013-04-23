function [ m ] = setbonds( m, blength )
%SETBONDS [ m ] = setbonds( m, blength )
% where m is a mat struct
    [ b, naa, pidx ] = mat2bond( m );

    bnorm = sqrt(sum(b.*b,2));

    bn = b./repmat(bnorm,1,3);
    bnew = bn*blength;
    shift = bnew - b;

    for i=2:naa,
        for j=pidx{i},
            m.atomr(j).x = m.atomr(j).x + sum(shift(1:i-1,1));
            m.atomr(j).y = m.atomr(j).y + sum(shift(1:i-1,2));
            m.atomr(j).z = m.atomr(j).z + sum(shift(1:i-1,3));
        end
    end


end
