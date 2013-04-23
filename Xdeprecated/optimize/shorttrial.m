% trying ampl program with all short proteins
% there are several weird errors (like 1 in 5 of these proteins have
% errors) 


load shortps % loads trial - indices of all short proteins (bw 5 and 25 aa's)
s = size(trial,2);
GDTs = zeros(1,s);
%newclasses = 1;

% parameters for runopt: (runopt calls ampl)
probname = 'tester'; algnum = 1; maxit = 3000; clco = 1; hp = 1; initials = 0;
avginitials = 1; thfrags = 1; zlu = 1; displaymode = 0; individualopt = 0;

for i = 4:49
    i
    p = trial(i)
    a = data{p}.seq;
    h = find(a > 20)
    if isempty(h) == 0
        continue;
    end
    size(a)
    GDTs(i) = runopt(data,class_info,newclasses,confus,probname,p,algnum,maxit,clco,hp,initials,avginitials,thfrags,zlu,displaymode,individualopt);
end