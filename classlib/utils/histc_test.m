%x = randn(1000,1);
%y = randn(1000,1)+1;

ed = -1:0.2:2;

n = histc([x y],ed);

subplot(2,2,1)
bar(ed,n,'stacked');
subplot(2,2,2)
bar(ed,n(:,1),'histc');

subplot(2,2,3)
bar(ed+0.1,n(:,1),1);
xlim([ed(1)-0.1 ed(end)-0.1]);
subplot(2,2,4)
bar(ed+0.1,n(:,2),1);
xlim([ed(1) ed(end)]);

figure;
bar_histc(ed,n,1);