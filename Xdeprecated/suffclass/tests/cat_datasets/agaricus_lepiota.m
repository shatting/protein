function [ data ] = agaricus_lepiota()

load agaricus-lepiota-cell;

mat = zeros(length(mushroomscell),23);

for i=1:length(mushroomscell),
    cur = mushroomscell(i);
    cll = regexp(cur, ',', 'split');
    data(i,:)=int32(char(cll{1}))';
end

for c = 1:23,
    col = data(:,c);
    [u, m ,n] = unique(col);
    data(:,c) = n;
end

% swap class to last place
data(:,[1,end]) = data(:,[end,1]);

end
