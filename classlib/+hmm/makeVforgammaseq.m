function [Vr,Vl] = makeVforgammaseq(gamma,confus,newclasses)
% seq = data{p}.seq 
% makes the potential 'matrix' V used in dynopt for a protein sequence seq
% Vr is a 3*9*(n-3) matrix used in dynoptR, where si = ci + 3xi -3 and
%                                                 xi = ci+1 + 3ti - 3 (xi in
%                                                 1..9)
% Vl is a 3*21*(n-3) matrix used in dynoptL, where si = xi + 3ci+1 - 3 and
%                                                  xi = ci + 9ti - 9 
%T and S are the tables for dynoptR (T) and dynoptL (S)



n = length(gamma); %technically n - 3, but whatever

if newclasses == 0
    Vr = zeros(3,9,n); %the potential matrix for dynoptR
    Vl = zeros(3,9,n); %the potential matrix for dynoptL

    xs = [1,2,3,10,11,12,19,20,21]; %the values of x for Vl


    for i = 1:n
        for u = 1:3
            for x = 1:9
                if gamma(i) ~= 0
                    if confus(gamma(i),u + 3*x - 3) ~= 0
                       Vr(u,x,i) = -log(confus(gamma(i),u + 3*x - 3));
                  %    Vl(u,x,i) = -log(confus(gamma(i),u + 3*x - 3));
                  %  else
                  %      disp('yuck - confusion matrix has a 0 here!');
                  %      p,u,x
                    end
                    %potential is V(gamma,sec) = -log(confus(gamma,s), s = c + 3x - 3 
                    %Vl(u,x,i) = -log(confus(gamma(i),x + 3*u - 3));
                    %be careful! the sequence goes in the reversed direction for
                    %gamma
                    if confus(gamma(i),3*u + xs(x) - 3) ~= 0
                        Vl(u,x,i) = -log(confus(gamma(i),3*u + xs(x) - 3));
                   % else
                   %     disp('yuck - confusion matrix has a 0 here!');
                   %     p,u,x
                    end
                else % if gamma(i) = 0 %this is a very big decision
                    Vl(u,x,i) = 0;
                    Vr(u,x,i) = 0; %setting V to zero causes the potential to be higher, making this sequence less likely
                end
            end
        end
    end

else %newclasses = 1
    Vr = zeros(1,4,n);
    Vl = zeros(1,4,n);
    
    for i = 1:n
        for u = 1:4
            if gamma(i) ~= 0
                    if confus(gamma(i),u) ~= 0
                       Vr(1,u,i) = -log(confus(gamma(i),u));
                       Vl(1,u,i) = -log(confus(gamma(i),u));
                    end
            end   
        end
    end
end %end newclasses
            
