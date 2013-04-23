%if ~ask('do optimization tests?'), return, end;

p = 323;
dprintf('Testing prediction on protein %i from database.',p);

if 1, %proteinnumber, database, new_deal_data, problemname, confus, potential, sbounds
    options = struct;
    options.HMM_numseqs = 19;
    optipot = optipotential_load('hel27_ndgold');
    [rmses,bests,bestpre] = predict(p,optipot,'hel27_ndgold',options);
end

if 0,
    disp('run optimization once, to see if everything is working!');
    runopt(data,class_info,confus,'tester',p,1,1000,0)
end
% trial of ampl program with all short p's:
if 0,
    shorttrial
    disp('sometimes shorttrial has errors- I think that they are due to proteins');
    disp('with aa > 20 or confus = 0 or other such problems');
end