g = whos('global');
disp('------- global strings ----------');
for i=1:length(g),
    if (strcmp(g(i).class,'char'))
        disp(sprintf('global %s = %s',g(i).name,eval(g(i).name)));
    end
end

clear g; clear i;