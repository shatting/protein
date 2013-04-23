function [S,SV,CL,jsofeachs,endj,endi,pastis,CLcol,jeq1,L] = getnbest(gammaseq,sgconfus,nfound,CL,jsofeachs,endj,endi,past,CLcol,S,SV,jeq1,L)
% this function should provide the same result as sequencesforp, but should
% work by finding exactly one sequence at a time, starting with the very
% best. (sequencesforp finds the nfound best sequences simultaneously). 
% the benefit of this method is that it stores all the info. calculated in 
% order to find each sequence, so that if more sequences are desired, it
% does not repeat any calculations. it is based on an article: "Sequentially
% Finding the N-Best List in HMMS"
% so, this will save the potentials of tons of sequences in a matrix called
% C. really, each entry of C just saves the best potential of all sequences
% in the set B(i) ~ {s|s_{1} = sb_{1}, s_{2}=sb_{2},...s_{i-1}=sb_{i-1}, s_{i}
% ~= sb_{i}}, where sb is one of the sequences already determined to be the
% nth best. (nth = first, second, third,..)
% then the minimum of C is found, and the sequence corresponding to that
% potential is derived, partially through its relationship to the sequence
% sb, and partially with the values Rk made by dynopt. 

%nums = nums;

if nargin < 8
    newsequences = 0;
else
    newsequences = 1;
    pastis = past;
end

% only hefg for now
newclasses = 1;
nfrags = length(gammaseq);

nums = size(sgconfus,2);

%make some matrices to reduce calculation time
% store -log(optidata.sgconfus) for less calculations.
minuslogconfus = -log(sgconfus);

% make the feasible pair matrix - should save a lot of time as opposed to
% using isfeasible
if newclasses == 0
    feasiblepairmatrix = makefeasiblepairmatrix;
else 
    feasiblepairmatrix = ones(nums); % all transitions allowed
end

%make a sys10toX matrix, since that is called so often
if newclasses == 0
    tentoXconverter = zeros(nums,3); % columns contain c_i, cp_i, t_i for rows i = 1:nums
    for s = 1:nums
        a = sys10toX(s,[3 3 3]);
        tentoXconverter(s,:) = a;
    end
end


%start by finding the minimal sequence:
[Vmin,Vminf,Rk,Lj,mins] = hmm.findpartialsumsforgammaseq(gammaseq,sgconfus,newclasses);
s1 = mins(:,1);
clear Lj;

if 1 % newclasses == 1
    Rkmin = min(Rk,[],1);
    Rk = Rkmin;
end

    
if newsequences == 0
    %make the matrix, where the nfound minimal sequences go
    S = zeros(nfound,nfrags);
    S(1,:) = s1;
    % the matrix, where the potentials of the minimal sequences go
    SV = zeros(nfound,1);
    SV(1) = Vmin;
    psums = zeros(nfound,3); %one column for psummin, one column for , one for sumb4
    psums(1,:) = [Vmin,Vmin,Vmin];

    %make the matrix CL, the 'candidate list', where the minimal potentials for
    %each Ai go
    CL = inf(nfrags,nfound);
    CLcol = 1; %CL column currently being developed


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Round 1: find s2 (the second best sequence!) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % now find first set of Ai's
    % or maybe just find their potentials
    % start with A(1), since it's a bit different
    % okay, so, in A(1), I'm supposed to store all possible sequences s with
    % s(1) ~= s1(1). But that's way too many sequences!!! Using Rk, I can
    % easily calculate the minimal sequence with s(1) ~= s1(1). So I'll do
    % that:
    % calculate the potential of s with s(1) ~= s1(1), for all s (except s1(1))
    A1 = inf(nums,1);
    sno = s1(1);
    if sno == 1
        for s = 2:nums
            if newclasses == 0
                cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                Rkindex = cpi; %minimum for s_(j+1).c = cpi 
                A1(s) = minuslogconfus(gammaseq(1),s) + Rk(Rkindex,2);
            else
                A1(s) = minuslogconfus(gammaseq(1),s) + Rkmin(2);
            end

        end
    elseif ((1 < sno) && (sno < nums))
        for s = 1:(sno-1)
            if newclasses == 0
                cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                Rkindex = cpi; %minimum for s_(j+1).c = cpi 
                A1(s) = minuslogconfus(gammaseq(1),s) + Rk(Rkindex,2);
            else
                A1(s) = minuslogconfus(gammaseq(1),s) + Rkmin(2);
            end
        end
        for s = (sno+1):nums
            if newclasses == 0
                cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                Rkindex = cpi; %minimum for s_(j+1).c = cpi 
                A1(s) = minuslogconfus(gammaseq(1),s) + Rk(Rkindex,2);
            else
                A1(s) = minuslogconfus(gammaseq(1),s) + Rkmin(2);
            end
        end
    else
        for s = 1:26
            if newclasses == 0
                cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                Rkindex = cpi; %minimum for s_(j+1).c = cpi 
                A1(s) = minuslogconfus(gammaseq(1),s) + Rk(Rkindex,2);
            else
                A1(s) = minuslogconfus(gammaseq(1),s) + Rkmin(2);
            end   
        end
    end % end if sno ==...
    CL(1,CLcol) = min(A1); %now add entry to C! 

    % so, now I have to do the same thing for i = 2:nfrags. it's trickier now, bc
    % the s_(1),...s_(i-1) = s1_(1),...s1_(i-1), and s_(i) ~+ s1_(i)

    % need to calculate potential of s1_(1)...s1_(i-1) in loop. this is
    % actually = Lk(i-1), since s1 = mins, but so that this first step is
    % compatible with all the other steps, I will add it up in the program. but
    % i can check if it is equal to Lk(i-1)! 
    sumb4 = 0;

    for i = 2:nfrags-1
        A = inf(nums,1);  
        % the new s_(i) has to be compatible with s1_(i-1)
        sminus = s1(i-1);
        % the new s_(i) can't be s1_(i)
        sno = s1(i);
        % need to know potential of s1_(1)...s1_(i-1)  
        sumb4 = sumb4 + minuslogconfus(gammaseq(i-1),s1(i-1));

        if sno == 1
            for s = 2:nums
                feas = feasiblepairmatrix(sminus,s);  %the new s_(i) has to be compatible with s1_(i-1)
                if feas == 1 % if sminu --> s possible
                    if newclasses == 0
                        cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                        Rkindex = cpi; %minimum for s_(j+1).c = cpi   
                       A(s) = sumb4 + minuslogconfus(gammaseq(i),s) + Rk(Rkindex,i+1);
                    else
                        A(s) = sumb4 + minuslogconfus(gammaseq(i),s) + Rkmin(i+1);
                    end
                end
            end
        elseif ((1 < sno) && (sno < nums))
            for s = 1:(sno-1)
                feas = feasiblepairmatrix(sminus,s);  %the new s_(i) has to be compatible with s1_(i-1)
                if feas == 1 % if sminu --> s possible
                    if newclasses == 0
                        cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                        Rkindex = cpi; %minimum for s_(j+1).c = cpi 
                        A(s) = sumb4 + minuslogconfus(gammaseq(i),s) + Rk(Rkindex,i+1);
                    else
                        A(s) = sumb4 + minuslogconfus(gammaseq(i),s) + Rkmin(i+1);
                    end
                end
            end
            for s = (sno+1):nums
                feas = feasiblepairmatrix(sminus,s);  %the new s_(i) has to be compatible with s1_(i-1)
                if feas == 1 % if sminu --> s possible
                    if newclasses == 0
                        cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                        Rkindex = cpi; %minimum for s_(j+1).c = cpi 
                        A(s) = sumb4 + minuslogconfus(gammaseq(i),s) + Rk(Rkindex,i+1);
                    else
                        A(s) = sumb4 + minuslogconfus(gammaseq(i),s) + Rkmin(i+1);
                    end
                end
            end
        else
            for s = 1:size(feasiblepairmatrix,2)
                feas = feasiblepairmatrix(sminus,s);  %the new s_(i) has to be compatible with s1_(i-1)
                if feas == 1 % if sminu --> s possible
                    if newclasses == 0
                        cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                        Rkindex = cpi; %minimum for s_(j+1).c = cpi 
                        A(s) = sumb4 + minuslogconfus(gammaseq(i),s) + Rk(Rkindex,i+1);
                    else
                        A(s) = sumb4 + minuslogconfus(gammaseq(i),s) + Rkmin(i+1);
                    end
                end
            end
        end % end if sno ==...
    CL(i,CLcol) = min(A);

    end %end for i = 2:nfrags-1

    % now: i = nfrags
    i = nfrags;
    A = inf(nums,1);  
    sumb4 = sumb4 + minuslogconfus(gammaseq(i-1),s1(i-1));
        % the new s_(i) has to be compatible with s1_(i-1)
        sminus = s1(i-1);
        % the new s_(i) can't be s1_(i)
        sno = s1(i);

        if sno == 1
            for s = 2:nums
                feas = feasiblepairmatrix(sminus,s);  %the new s_(i) has to be compatible with s1_(i-1)
                if feas == 1 % if sminu --> s possible
                       A(s) = sumb4 + minuslogconfus(gammaseq(i),s);
                end
            end
        elseif ((1 < sno) && (sno < nums))
            for s = 1:(sno-1)
                feas = feasiblepairmatrix(sminus,s);  %the new s_(i) has to be compatible with s1_(i-1)
                if feas == 1 % if sminu --> s possible
                    A(s) = sumb4 + minuslogconfus(gammaseq(i),s);
                end
            end
            for s = (sno+1):nums
                feas = feasiblepairmatrix(sminus,s);  %the new s_(i) has to be compatible with s1_(i-1)
                if feas == 1 % if sminu --> s possible
                    A(s) = sumb4 + minuslogconfus(gammaseq(i),s);
                end
            end
        else
            for s = 1:26
                feas = feasiblepairmatrix(sminus,s);  %the new s_(i) has to be compatible with s1_(i-1)
                if feas == 1 % if sminu --> s possible
                    A(s) = sumb4 + minuslogconfus(gammaseq(i),s);
                end
            end
        end % end if sno ==...
    CL(i,CLcol) = min(A);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % so now that I've found all the potentials, I choose the minimal %
    % potential and find the corresponding sequence -  s2!            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %[pot2,ind] = min(CL(:,CLcol));
    [pot2,ind,j] = minindex(CL(:,CLcol));
    % the index tells us up to which index s2 is identical to s1 (actually sj) - nifty!
    % s_(ind) is the first index different than s1
    % calculate potential while finding sequence
    % I think the potential should actually be saved somewhere! (really it's
    % Lk(i), but when I'm dealing with different iterations (it's only Lk for step 1-2), it won't work that
    % way, so I want to do it the hard way now, and then economize later)
    sumb4 = 0;
    for i = 1:ind-1
        S(2,i) = S(1,i);
        sumb4 =  sumb4 + minuslogconfus(gammaseq(i),S(2,i));
    end
    % now need to find the end of the sequence! a bit more work
    % sequence has to be feasible with minimal potential
    % first s2_(ind)

    psummin = inf;

    for s = 1:nums
        if s == S(1,ind)
            continue;
        end
        if ind > 1
          feas = feasiblepairmatrix(S(2,ind-1),s);
        else
           feas = 1; % every s is feasible at the first position!
        end
        if feas == 1
            if newclasses == 0
                cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                Rkindex = cpi; %minimum for s_(j+1).c = cpi 
            else
                Rkindex = 1;
            end
            if ind < nfrags
               psum = sumb4 + minuslogconfus(gammaseq(ind),s) + Rk(Rkindex,ind+1);
            else
                psum = sumb4 + minuslogconfus(gammaseq(ind),s);
            end
            if psum < psummin  
                psummin = psum;
                sind = s;
            end
        end
    end
    S(2,ind) = sind;
    % check! 
    %psummin, pot2
    sumb4 =  sumb4 + minuslogconfus(gammaseq(ind),S(2,ind));


    % now the rest of the sequence!
    for i = ind+1:nfrags
       psummin = inf;

       for s = 1:nums
          feas = feasiblepairmatrix(S(2,i-1),s);
          if feas == 1
              if newclasses == 0
                  cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                  Rkindex = cpi; %minimum for s_(j+1).c = cpi 
              else
                  Rkindex = 1;
              end
            if i < nfrags
               psum = sumb4 + minuslogconfus(gammaseq(i),s) + Rk(Rkindex,i+1);
            else
                psum = sumb4 + minuslogconfus(gammaseq(i),s);
            end
            if psum < psummin  
                psummin = psum;
                si = s;
            end
          end
       end

    S(2,i) = si;
    %psummin, pot2
    sumb4 =  sumb4 + minuslogconfus(gammaseq(i),S(2,i));

    end

    SV(2) = sumb4;
    psums(2,:) = [psummin,pot2,sumb4];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Round 2 -end: find s3 - s_nfound - the rest of the sequences! %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % tricky, tricky
    % the indices don't match the earlier indices
    i = j;
    j = ind;

    %S now has 2 rows, so L = 2, and we are looking for s^{L = 3}

    %usedsjs = cell(nfrags,nfound); % to keep track of sj's used coming from a certain position in CL
    %usedsjs{j,i } = S(1,j);
    %usedsjs{j,i} = [usedsjs{j,i},S(2,j)];
    % make a matrix to keep track of which CL(.,i) each s^{L} was in
    jsofeachs = zeros(nfound - 1,1);
    jsofeachs(1) = j;  %the i belonging to s^{2} - so s^{2}_{i} ~= s^{all other s's with this i}_{i}
    jeq1 = [S(1,1)];  % keep track of all s's in A(i,1) - the next time s in in A(i,1), it can't have s_(1) = any of these
    pastis = zeros(nfound,1);    % keep track of the column i that each sequence was in. the future s values depend on where the s with j's before this s were, and where the ones be4 them were
    pastis(2,1) = i;

    a = 3; % index where we start loop ( a = number of next sequence we are looking for)
else % end if newsequences = 0
    a = L+1;
    CL = [CL,inf(nfrags,nfound)];
    nfound = a + nfound - 1;
    j = endj;
    i = endi;
end

for L = a:nfound  %the index of the sequence we are looking for: s^{L}
    
    CLcol = CLcol + 1; % we are now working on the next column of CL!!!
    % start by figuring out which partition A of the set of sequences contains
    % the sequence just found
    % we did [pot2,ind,j] = minindex(CL(:,CLcol)), so this is equivalent to what will be
    % done at the end of each of these rounds: [potn,j,i] = minindex(CL)
    % so we know that the last sequence came from the jth row of CL(:,i), 
    % so we now have to break up the set of sequences D whose minimal
    % potential is CL(j,i)
    %s^{L-1} in C(j,i) --> s^{L-1}_{1} = s^{i}_{1},s^{L-1}_{2} =
    %s^{i}_{2},..s^{L-1}_{j-1} = s^{i}_{j-1}, and s^{L-1}_{j} ~= s^{i}_{j}
    
    % need to find the minimal potential of D_{r} for r = j,...nfrags (D_{r} =
    % inf, for r = 1,..,j-1
    
    % potential of already fixed s_{i}'s - this should actually be stored
    % somewhere!
    sumb4 = 0;
    for k = 1:j-1
        sumb4 = sumb4 + minuslogconfus(gammaseq(k),S(L-1,k));
    end
    sumfirst = sumb4;
    
    %important!! set C(j,i) = inf!!
    CL(j,i) = inf;
    
    % start with D_{j}, definied slightly differently
    % D_{j} = s, up until s_{j}, s_{k} = s(L-1)_{k}
     % check which s's have been used in previous sequences with same
     % beginning part
     if j == 1
         jeq1 = [jeq1,S(L-1,1)];
         sno = jeq1; % s_(1) can't be any of s's with j = 1, and can't be S(1,1)
     else % j > 1
        sinos = find(jsofeachs == j);  %((s == S(L-1,j)) || (s == S(i,j))) 
        sno = zeros(1,size(sinos,2));
        for r = 1:size(sinos,1)
           sno(r) = S(sinos(r)+1,j); %have to look at sequences found with these j's
        end
        sno = [sno,S(L-1,j),S(i,j)];
     end
     
     % add sno from past i's. call this snow for now, bc I think it's going to
     % overlap, and then I can improve sno to that there are no repeated
     % indices
     r = i;%L-1;
     snow = S(r,j);
     while r> 1
         r = pastis(r);
         probs = 0; %make sure the sequences are the same up until j
         for b = 1:j-1
            if S(r,b) ~= S(L-1,b)
                probs = probs + 1;
            end
         end
         if probs == 0
             snow = [snow,S(r,j)];
         end
     end
     snow = [snow,S(L-1,j),S(i,j)];
      if j == 1
         snow = [snow,sno];
     end
     
    %minVDj = inf;
    VDjs = inf(1,nums);
    for s = 1:nums
        % avoid s's already used in sequences starting this way
%         if isempty(find(sno == s)) ~= 1
%              continue;
%         end
        if isempty(find(snow == s)) ~= 1
             continue;
        end
        % s^{L}_{j} also can't be s^{L-1)_{j}
       % if ((s == S(L-1,j)) || (s == S(i,j)))
       %     continue;
       % end
        % check if s is feasible
        if j > 1
            feas = feasiblepairmatrix(S(L-1,j-1),s);
        else
            feas = 1;
        end
        if feas == 1
            if newclasses == 0
                cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                Rkindex = cpi; %minimum for s_(j+1).c = cpi 
            else
                Rkindex = 1;
            end
            if j < nfrags
                VDjs(s) = sumb4 + minuslogconfus(gammaseq(j),s) + Rk(Rkindex,j+1);
            else 
                VDjs(s) = sumb4 + minuslogconfus(gammaseq(j),s);
            end
%             if VDj < minVDj
%                 minVDj = VDj;
%             end
        end %if feas == 1
    end %for s = 1:nums
    CL(j,CLcol) = min(VDjs);
    
    for k = j+1:nfrags
        % look at each sequence in D, and see if pot smaller than min!
        %minVDk = inf;
        Vdks = inf(1,nums);
        % increase the sum of fixed sequence members:
        sumb4 = sumb4 + minuslogconfus(gammaseq(k-1),S(L-1,k-1));
        
        for s = 1:nums
             % check if s is feasible and not s^(L-1)_(k)
            if s == S(L-1,k)
                continue;
            end
            feas = feasiblepairmatrix(S(L-1,k-1),s);
            if feas == 1
                if newclasses == 0
                    cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                    Rkindex = cpi; %minimum for s_(j+1).c = cpi  
                else
                    Rkindex = 1;
                end
                if k < nfrags
                    Vdks(s) = sumb4 + minuslogconfus(gammaseq(k),s) + Rk(Rkindex,k+1);
                else
                    Vdks(s) = sumb4 + minuslogconfus(gammaseq(k),s);
                end
               % if VDk < minVDk
               %     minVDk = VDk;
               % end
            end %if feas == 1
        end %for s = 1:nums
        
        CL(k,CLcol) = min(Vdks);
     end %for k = j+1:nfrags
         
     
     
     
     %% now find minimal position in CL -  the potential of the next
     %% minimal sequence, s^{L}
   %  jold = j;
   %  iold = i;
     [potn,j,i] = minindex(CL);
     
     pastis(L,1) = i; %keep track of which i this L was in
     % also need to know the i's of other L's with this j
     %for u = 2:L-1
     %     if jsofeachs(u-1) == j
     %        pastis(L,u) = u;
     %    end
     %end
     
     % the same minimum shows up several times in CL!!! is that okay? it
     % seems correct...
     % now use potn to find s^{L}: (actually, we don't ever use potn, just
     % its indices in CL ! (is that okay?)
     
     %uh-oh! it doesn't seem right that I'm doing this twice in each
     %iteration
     % check which s's have been used in previous sequences with same
     % beginning part
%      sinos = find(jsofeachs == j);  %((s == S(L-1,j)) || (s == S(i,j))) 
%      sno = zeros(size(sinos));
%      for r = 1:size(sinos,1)
%          sno(r) = S(sinos(r),j);
%      end
     
     
     if j > 1
       S(L,1:j-1) = S(i,1:j-1); %the first part of s^{L} is identical to s^{i}
     end
     sumb4 = 0; %not sumfirst from before, that was a different j; % the potential of the first part of the sequence
     % now have to find the rest of the terms!
     for k = 1:j-1
         sumb4 = sumb4 + minuslogconfus(gammaseq(k),S(L,k));
     end
     
     if j == 1
         sno = [jeq1,S(1,j)]; % s_(1) can't be any of s's with j = 1, and can't be S(1,1)
     else % j > 1
         %then s_{j} can't be s^{L-1}_(j) or s^{i)_{j}
         sinos = find(jsofeachs == j);  %((s == S(L-1,j)) || (s == S(i,j))) 
        sno = zeros(1,size(sinos,2));
        for r = 1:size(sinos,1)
           sno(r) = S(sinos(r)+1,j); %have to look at sequences found with these j's
        end
        sno = [sno,S(i,j)]; %do I need S(L-1,j), too?
     end
     
     r = i;%L-1;
     snow = S(r,j);
    while r> 1
         r = pastis(r);
         probs = 0; %make sure the sequences are the same up until j
         for b = 1:j-1
            if S(r,b) ~= S(L,b)
                probs = probs + 1;
            end
         end
         if probs == 0
             snow = [snow,S(r,j)];
         end
     end
     snow = [snow,S(i,j)];
     if j == 1
         snow = [snow,sno];
     end
     
     % start with s_{j} (can't be equal to S^{i}_{j} or S^{L-1}_{j}
     % it also can't be equal to S^{ }_{j} for any other sequences that
     % came from C(j,i)!!!! so need to make a list of sj's for C(j,i)
     psummin = inf;
     for s = 1:nums
%         if isempty(find(sno == s)) ~= 1
%              continue;
%         end
        if isempty(find(snow == s)) ~= 1
             continue;
        end
        % s^{L}_{j} also can't be s^{L-1)_{j}
%         if ((s == S(L-1,jold)) || (s == S(iold,jold)))
%             continue;
%         end
        if j > 1
          feas = feasiblepairmatrix(S(L,j-1),s);
        else
           feas = 1; % every s is feasible at the first position!
        end
        if feas == 1
            if newclasses == 0
                cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                Rkindex = cpi; %minimum for s_(j+1).c = cpi 
            else
                Rkindex = 1;
            end
            if j < nfrags
               psum = sumb4 + minuslogconfus(gammaseq(j),s) + Rk(Rkindex,j+1);
            else
                psum = sumb4 + minuslogconfus(gammaseq(j),s);
            end
           % if psum < psummin  % try next: if psum == potn
           %     psummin = psum;
           %     sj = s;
           % end
           if (psum - potn) < 10^(-10)
               sj = s;
               continue;
           end
        end
     end %end for s 
    S(L,j) = sj;
    % update list of sj's coming from C(j,i) - make a cell of size (nfrags,nfound),
    % and in each entry have a list of sj's used!
    jsofeachs(L-1) = j;
    
    sumb4 =  sumb4 + minuslogconfus(gammaseq(j),S(L,j));
    
    % check!
    
    %now find the rest of the s_{k}'s (k = j+1 : nfrags)
    for k = j+1:nfrags
       psummin = inf;

       for s = 1:nums
          feas = feasiblepairmatrix(S(L,k-1),s);
          if feas == 1
              if newclasses == 0
                  cpi = tentoXconverter(s,2); %a(1); cpi = a(2); 
                  Rkindex = cpi; %minimum for s_(j+1).c = cpi 
              else
                  Rkindex = 1;
              end
            if k < nfrags
               psum = sumb4 + minuslogconfus(gammaseq(k),s) + Rk(Rkindex,k+1);
            else
                psum = sumb4 + minuslogconfus(gammaseq(k),s);
            end
            %if psum < psummin  
            %    psummin = psum;
            if (psum - potn) < 10^(-10)
                sk = s;
            end
          end
       end

    S(L,k) = sk;
    sumb4 =  sumb4 + minuslogconfus(gammaseq(k),S(L,k));

    end %end for k = j+1:nfrags
    
    SV(L) = sumb4; 
    psums(L,:) = [psummin,potn,sumb4];
    
    
end %end for L = 3:nfound

endj = j;
endi = i;

return
%% CHECK!!!
A = zeros(nfound,nfound);
for k = 1:nfound
  for l = 1:nfound
    A(k,l) = opti.percentsequenceidentity(S(k,:),S(l,:));
    if A(k,l) > .999999
        if k ~= l
          k
          l
          A(k,l)
        end
    end
  end
end





