function feasiblepairmatrix = makefeasiblepairmatrix
%matrix which tells if it is possible to go from s1 (rows) to s2 (columns)
% 1 = possible, 2 = not possible

nums = 27;
feasiblepairmatrix = zeros(nums,nums);

for i = 1:nums
    for j = 1:nums
        problems = issequencefeasible([i j]);
        if problems == 0 %pair is possible!
           feasiblepairmatrix(i,j) = 1;
        else
            feasiblepairmatrix(i,j) = 0;
        end
    end
end