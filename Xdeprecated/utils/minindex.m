function [mina,i,j] = minindex(a)
% finds the minimum value of an (m x n) matrix, plus the row and column number
% so that a(i,j) = mina = min(min(a))
mina = min(min(a));
[i,j] = find(a == mina);
i = i(1);
j = j(1);


