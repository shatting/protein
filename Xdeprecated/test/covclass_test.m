
numpoints = 100;

p1 = rand(numpoints,2);
p2 = rand(numpoints,2)+2;

%p1(1,:) = [0,0];
%p2(numpoints,:) = [2,2];

p=[p1;p2];


[pot,cl,confus] = covclass([p1;p2],[1;zeros(numpoints-1,1);zeros(numpoints-1,1);2])

%covclass_test_plot(1,p,cl,pot);
