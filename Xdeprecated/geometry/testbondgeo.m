
%if 0;

    maxb = -inf;
    problems = [];
    for p = 1:1334
        if mod(p,100) == 1
            p
            maxb
        end

        bond = data{p}.bond;

        bond = rotatematrix(bond); %normalizes and rotates bond


        geo = bond2geo(bond);

        bondnew = geo2bond(geo);


        maxbond = max(max((bondnew - bond)));
        if maxbond > maxb
            maxb = maxbond;
        end
        if maxbond > 10^(-5)
            problems = [problems,p];
        end

    end

%end %end if 0

% this part tests to see if tpo, to can be calculated from tqu, based on
% tcl
if 0
problems = [];
for p = 1:1334
    if mod(p,100) == 0
        p
        problems
    end

    bond = data{p}.bond;
    bond = rotatematrix(bond); %normalizes and rotates bond
    geo = bond2geo(bond);

    A = zeros(size(geo.co,1),6);
    for i = 1:size(geo.co,1)
        q = geo.tqu(i);
        t0 = geo.t0(i);
        if geo.tcl(i) == 1
            % in this case, tqu = tsin = to
            t = q/t0;
            tp = sqrt(1 - (q/t0)^2);
            %[t tp t*geo.t0(i) tp*geo.t0(i) q q*geo.t0(i) geo.to(i) geo.tpo(i) geo.to(i)*geo.t0(i) geo.tpo(i)*geo.t0(i)]
        elseif  geo.tcl(i) == 2
            % in this case, tqu = tcos = tpo
             tp = q/t0;
             t = sqrt(1 - (q/t0)^2);
            %[t tp t*geo.t0(i) tp*geo.t0(i) q q*geo.t0(i) geo.to(i) geo.tpo(i) geo.to(i)*geo.t0(i) geo.tpo(i)*geo.t0(i)]

        elseif  geo.tcl(i) == 3
            % in this case, tqu = tsin = to 
             t = q/t0;
             tp = -sqrt(1 - (q/t0)^2);

            %[t tp t*geo.t0(i) tp*geo.t0(i) q q*geo.t0(i) geo.to(i) geo.tpo(i) geo.to(i)*geo.t0(i) geo.tpo(i)*geo.t0(i)]
        elseif  geo.tcl(i) == 4
            % in this case, tqu = tcos = tpo
             tp = q/t0;
             t = -sqrt(1 - (q/t0)^2);

            %[t tp t*geo.t0(i) tp*geo.t0(i) q q*geo.t0(i) geo.to(i) geo.tpo(i) geo.to(i)*geo.t0(i) geo.tpo(i)*geo.t0(i)]
        end
        A(i,:) = [geo.to(i),t,geo.tpo(i),tp,geo.tcl(i),q];
    end
    if (max(A(:,1) - A(:,2)) > 10^-5) || (max(A(:,3) - A(:,4)) > 10^-5)
        problems = [problems,p];
    end

end %end for p

end % end if 0 

