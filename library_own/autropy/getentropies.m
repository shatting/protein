function [ entout, delta ] = getentropies( Vgd )
%GETENTROPIES Summary of this function goes here
%   Detailed explanation goes here

setup_cvx


ncl = size(Vgd,1);

% empty classes. cvx doesnt like inf in constraints
idxninf = max(triu(Vgd,1),[],2) ~= inf;
nclninf = sum(idxninf);

Vgd = Vgd(idxninf,idxninf);

cvx_begin
    variables delta ent(nclninf);
    maximize(delta);
    subject to
        for g=1:nclninf,
            for d=1:nclninf,
             if d~=g,
                ent(g) - ent(d) + delta <= Vgd(g,d);  
             end;
            end
        end
        ent(nclninf) == 0;
cvx_end

% only one loop try
% triag = triu(ones(nclninf),1);
% lin = zeros(nclninf);
% lin(triag==1) = 1:(nclninf)*(nclninf-1)/2;
% 
% S = sparse(lin);
% [g,d,s] = find(S);
% 
% nineq = length(g);
% 
% cvx_begin
%     variables delta ent(nclninf);
%     maximize delta;
%     subject to
%         for i=1:nineq,
%             ent(g(i)) - ent(d(i)) + delta <= Vgd(g(i),d(i));
%         end
%         ent <= 10;
%         ent(1) == 0;
% cvx_end
%----

entout = zeros(ncl,1) + inf;
entout(idxninf) = ent;

end
