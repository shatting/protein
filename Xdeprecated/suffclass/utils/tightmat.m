function tightmat( m, format, integerformat )
%TIGHTMAT print integer matrix in tight space

[a,nc] = size(m);


if nargin<2, 
    x = max(m,[],1)';
    y = ceil(log(x)/log(10));
    for i=1:nc,
        if (abs(y(i))<inf)
           s =  ['%' int2str(y(i)+2) '.f'];
        else
           s = ['%4.f'];
        end
        format{i} = s;
    end
end


s='';
for c=1:nc,
    if (iscell(format)),
        s=[s format{c}]; 
    elseif (max(abs(rem(m(:,c),1))) < 1e-10)
        s=[s integerformat];
    else
        s=[s format]; 
    end
end
s = [s '\n'];

dprintf(s,m');

end
