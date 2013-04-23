clf;

data = rand(5,5,10);

for i=1:size(data,1),
    i
    for j=1:size(data,2),
        histo(data(i,j,:),[i,j], sprintf('%i,%i,1',i,j));
    end;
end;
