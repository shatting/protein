% PLOT_HISTO_TEST Test plot_histo.
clf;
data = rand(5,5,10);
hold on;

for i=1:size(data,1),
    i
    for j=1:size(data,2),
        plot_histo(data(i,j,:),[i,j], sprintf('%i,%i,1',i,j));

    end;
end;
