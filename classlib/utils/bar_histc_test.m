% BAR_HISTC_TEST Test bar_histc.m

x = randn(1000,1);
y = randn(1000,1)+1; % shift by 1

edges = -1:0.2:2;

histdata = histc([x y],edges);

subplot(3,2,1)
bar(edges,histdata(:,1),'histc');
xlabel('data 1');

subplot(3,2,2)
bar(edges+0.1,histdata(:,2),1);
xlabel('data2');

subplot(3,2,[3,4])
bar(edges,histdata,'stacked');
xlabel('stacked');


subplot(3,2,[5,6]);
bar_histc(edges,histdata,0.9);
xlabel('bar\_histc: the bars are at the correct positions');