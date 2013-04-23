function [ C, Ceq ] = fun_nonlincon( x, sseq, sbounds )
% C(x) <= 0, Ceq(x) == 0

nfrags = length(sseq);

for f = 1:nfrags,
    ind = 3*(f-1)+1:3*(f-1) + 3 ;
    % f = 1 -> 1:3
    % f = 2 -> 4:6    
    x_i = x(ind)';
    
end

fval = sum(individualpots);

end
