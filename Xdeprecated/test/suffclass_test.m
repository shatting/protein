addpath '../suffclass';

tic
numpoints = 100;

p1 = rand(numpoints,2);
p2 = rand(numpoints,2)+2;

%p1(1,:) = [0,0];
%p2(numpoints,:) = [2,2];

p = [p1;p2];

suffs = cell(numpoints*2,1);

emptystat = suffstat(1,zeros(2,1));

for i=1:size(p,1),
    suffs{i} = suffstat(emptystat,p(i,:),1,1);
end


[pot,cl,confus] = suffcov(suffs,[1;zeros(numpoints-1,1);zeros(numpoints-1,1);2])

%covclass_test_plot(1,p,cl,pot);
