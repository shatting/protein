function problems = issequencefeasible(seq)
% checks if a class sequence is feasible (according to c_i, c_i+1)
% nums is the number of s classes (was 27)
% seq is a class sequence displayed horizontally
% checked this program, and all actual class sequences are feasible (good!)
% be very careful!! sys10toX only works for doubles!!!
% ( problems = 0 --> sequence is feasible)

sizes = length(seq);
 
    
%alternative, possibly speedier method 
problems = 0;
geo = sys10toX(seq,[3 3 3]);
for i = 2:sizes
    problems = problems + (geo(i-1,2) ~= geo(i,1));
end


%old, not so speedy method
    % if sizes == 2
%     a = seq(1);
%     b = seq(2);
% end
% cp = zeros(sizes,1);
% 
% if sizes == 2
%     geo = sys10toX(a,[3 3 3]);
% else
%     geo = sys10toX(seq(1),[3 3 3]);
% end
% 
% cp(1) = geo(2);
% problems = 0;
% 
% for i = 2:sizes
%     if sizes == 2
%        geo = sys10toX(b,[3 3 3]);
%     else
%         geo = sys10toX(seq(i),[3 3 3]);
%     end
%     cp(i) = geo(2);
%     if geo(1) ~= cp(i-1) %if c_i ~= cp_i-1
%         problems = problems + 1;
%     end
% end


