function starmat( m, fn )
%STARMAT Summary of this function goes here
%   Detailed explanation goes here

if nargin>1,
    m = fn(m);
end

[r,c] = size(m);

%t = repmat('X',1,c+1);
numbs = '';
for i=1:c, numbs = [numbs, sprintf('%i',mod(i,10))]; end;
    
t = [numbs, sprintf('\n')];
for i=1:r,
    for j=1:c,
       if m(i,j), s='*'; else s=' '; end
       t = [t, s];
    end
    t = [t, sprintf('%i\n',mod(i,10))];
end

t = [t, numbs];

disp(t);

end
