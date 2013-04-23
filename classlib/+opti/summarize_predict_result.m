function [ ] = summarize_predict_result( allpres, rmsds )
%VISUALIZE_PREDICTIONS Summary of this function goes here
%   Detailed explanation goes here

% 0: The final solution satisfies the termination conditions for verifying optimality.
% -100 to -199: A feasible approximate solution was found.
% -200 to -299: The code terminated at an infeasible point.
% -300: The problem was determined to be unbounded.
% -400 to -499: The code terminated because it reached a pre-defined limit.
% -500 to -599: The code terminated with an input error or some non-standard error.

its = zeros(length(allpres),1);
funcounts = its;
fvals = its;
exitflags = its;
times = its;

for i=1:length(allpres)
    its(i) = allpres(i).OUTPUT.iterations;
    funcounts(i) = allpres(i).OUTPUT.funcCount;
    fvals(i) = allpres(i).FVAL;
    exitflags(i) = allpres(i).EXITFLAG;
    %times(i) = allpres{i}.knitrotime;
end
dprintf('\n-- PREDICTION RUN RESULTS ---');
dprintf('%i predictions',length(allpres));
dminmax('iteration count',its);
dminmax('function count',funcounts);
dminmax('function value',fvals);
%dminmax('knitro time',times);
exitperfect = exitflags == 0;
exitok = exitflags <= -100 & exitflags >= -199;
exitunbounded = exitflags == -300;
exitinfeasible = exitflags <= -200 & exitflags >= -299;
exitknown = exitinfeasible | exitunbounded | exitok | exitperfect;

dprintf('exact solutions:     %3i',sum(exitperfect));
dprintf('approx solutions:    %3i',sum(exitok));
dprintf('unbounded problems:  %3i',sum(exitunbounded));
dprintf('infeasible problems: %3i',sum(exitinfeasible));
if (any(~exitknown))
    dprintf('other flags: \t\t%i',sum(~exitknown));
    disp(unique(exitflags(~exitknown)));
end

%dminmax('alternative rmsds',armsds);
dminmax('rmsd',rmsds);
dprintf('total function count: %i',sum(funcounts));

figure;
subplot(2,1,1);
%legend({pnormalized('rmsd','r',rmsds),pnormalized('itcount','b',its),pnormalized('fval','k',fvals)});%,pnormalized('ktime','g',times)});
legend({pmin('mrmsd','r',rmsds)});%,pnormalized('ktime','g',times)});
hold off;
subplot(2,1,2);
legend({pmin('fval','k',fvals)});
hold off;

end

function dminmax(item, vals)
    [m,mi] = min(vals);
    [mx,mxi] = max(vals);
    mn = mean(vals);
    sigm = std(vals);
    dprintf('%s in [%.2f (HMM %i), %.2f (HMM %i)] mean=%.2f +/- %.2f sigma',item,m,mi,mx,mxi,mn,sigm);            
end

function legend = pmin(item,style,vals)
    [~,maxidx] = min(vals);
    plot(vals,style);
    legend = sprintf('min %s = %.2f',item,vals(maxidx));
    hold on;
    plot(maxidx,vals(maxidx),[style,'o']);
end

function legend = pnormalized(item,style,vals)
    [~,maxidx] = max(abs(vals));
    plot(vals/vals(maxidx),style);
    legend = sprintf('1 ~ %s = %.2f',item,vals(maxidx));
    hold on;
end