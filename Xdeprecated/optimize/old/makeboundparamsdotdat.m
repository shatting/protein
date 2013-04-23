function makeboundparamsdotdat(name,nums,individualopt,clower,cupper,cplower,cpupper,zboth,zupper)
% if individualopt = 1, then make individual z bounds for all fragments
% still have to be programmed in:
% if individualopt = 2, then make individual c bounds for all fragments 
% if individualopt = 3, then make individual c and z bounds for all fragments
% for the individual bounds, the input should already be the new bounds

% sometimes we use zlower and zupper, sometimes just zboth
if nargin < 9
    zlu = 0;
else
    zlu = 1;
    zlower = zboth;
end

global potparamswzupperlower;

fid=fopen(name,'w'); 


    
%printf(fid,textv);
textv = ['#bound constraints for ci, ti dependent on si:\n'];

%%%%%%%%%%%%%%%%%%%%%%%%%

if ((individualopt == 2) || (individualopt == 3)) %individual c bounds
    sizep = size(clower,1);
    textv = [textv,'# individual bounds for c \n'];
    textv = [textv,'param seqlowerc := \n'];
    for i = 1:sizep - 1
        if (abs(clower(i)) ~= inf)
            textv = [textv,num2str(i),'  ',num2str(clower(i)),'\n'];
        else
            if clower(i) > 0
                textv = [textv,num2str(i),'  Infinity \n'];
            else
                textv = [textv,num2str(i),' -Infinity \n'];
            end
        end
    end
    if abs(clower(sizep)) ~= inf;
        textv = [textv,num2str(sizep),'   ',num2str(clower(sizep)),';\n'];
    else
        if clower(sizep) > 0
            textv = [textv,num2str(sizep),'  Infinity; \n'];
        else
            textv = [textv,num2str(sizep),' -Infinity; \n'];
        end
    end
    textv = [textv,'\n'];

    textv = [textv,'param sequpperc := \n'];
    for i = 1:sizep - 1
        if (abs(cupper(i)) ~= inf)
            textv = [textv,num2str(i),'  ',num2str(cupper(i)),'\n'];
        else
            if cupper(i) > 0
                textv = [textv,num2str(i),'  Infinity \n'];
            else
                textv = [textv,num2str(i),' -Infinity \n'];
            end
        end
    end
    if abs(cupper(sizep)) ~= inf;
        textv = [textv,num2str(sizep),'   ',num2str(cupper(sizep)),';\n'];
    else
        if cupper(sizep) > 0
            textv = [textv,num2str(sizep),'  Infinity; \n'];
        else
            textv = [textv,num2str(sizep),' -Infinity; \n'];
        end
    end
    textv = [textv,'\n'];,3

else %end individualopt == 2

%%%%%%%%%%%%%%%%%%%%%%%
% end of individual opt 2,3
%%%%%%%%%%%%%%%%%%%%%%%%

textv = [textv,'param sequpperc := \n'];
for i = 1:nums - 1
    if abs(cupper(i)) ~= inf;
        textv = [textv,num2str(i),'  ',num2str(cupper(i)),'\n'];
    else
        if cupper(i) > 0
            textv = [textv,num2str(i),'  Infinity \n'];
        else
            textv = [textv,num2str(i),' -Infinity \n'];
        end
    end
end
% i = nums:
if abs(cupper(nums)) ~= inf;
    textv = [textv,num2str(nums),'   ',num2str(cupper(nums)),';\n'];
else
    if cupper(nums) > 0
        textv = [textv,num2str(nums),'  Infinity; \n'];
    else
        textv = [textv,num2str(nums),' -Infinity; \n'];
    end
end
textv = [textv,'\n'];
textv = [textv,'param seqlowerc := \n'];
for i = 1:nums - 1
    if abs(clower(i)) ~= inf;
        textv = [textv,num2str(i),'  ',num2str(clower(i)),'\n'];
    else
        if clower(i) > 0
            textv = [textv,num2str(i),'  Infinity \n'];
        else
            textv = [textv,num2str(i),' -Infinity \n'];
        end
    end
end

% i = nums:
if abs(clower(nums)) ~= inf;
    textv = [textv,num2str(nums),'   ',num2str(clower(nums)),';\n'];
else
    if clower(nums) > 0
        textv = [textv,num2str(nums),'  Infinity; \n'];
    else
        textv = [textv,num2str(nums),' -Infinity; \n'];
    end
end
textv = [textv,'\n'];

end %%% end if individualopt == 2 or 3 (individual c bounds)

%%%%%%%%%%%%%%%%%%%
% end of c bounds, beginning of cp bounds
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%

if ((individualopt == 2) || (individualopt == 3)) %individual cp bounds
    sizep = size(cplower,1);
    textv = [textv,'# individual bounds for cp \n'];
    textv = [textv,'param seqlowercp := \n'];
    for i = 1:sizep - 1
        if (abs(cplower(i)) ~= inf)
            textv = [textv,num2str(i),'  ',num2str(cplower(i)),'\n'];
        else
            if cplower(i) > 0
                textv = [textv,num2str(i),'  Infinity \n'];
            else
                textv = [textv,num2str(i),' -Infinity \n'];
            end
        end
    end
    if abs(cplower(sizep)) ~= inf;
        textv = [textv,num2str(sizep),'   ',num2str(cplower(sizep)),';\n'];
    else
        if cplower(sizep) > 0
            textv = [textv,num2str(sizep),'  Infinity; \n'];
        else
            textv = [textv,num2str(sizep),' -Infinity; \n'];
        end
    end
    textv = [textv,'\n'];

    textv = [textv,'param sequppercp := \n'];
    for i = 1:sizep - 1
        if (abs(cpupper(i)) ~= inf)
            textv = [textv,num2str(i),'  ',num2str(cpupper(i)),'\n'];
        else
            if cpupper(i) > 0
                textv = [textv,num2str(i),'  Infinity \n'];
            else
                textv = [textv,num2str(i),' -Infinity \n'];
            end
        end
    end
    if abs(cpupper(sizep)) ~= inf;
        textv = [textv,num2str(sizep),'   ',num2str(cpupper(sizep)),';\n'];
    else
        if cpupper(sizep) > 0
            textv = [textv,num2str(sizep),'  Infinity; \n'];
        else
            textv = [textv,num2str(sizep),' -Infinity; \n'];
        end
    end
    textv = [textv,'\n'];,3

else %end individualopt == 2

    %%%%%%%%%%%%%%%%%%%%%%%
    % end of individual opt 2,3
    %%%%%%%%%%%%%%%%%%%%%%%%


    textv = [textv,'param sequppercp := \n'];
    for i = 1:nums - 1
        if abs(cpupper(i)) ~= inf;
            textv = [textv,num2str(i),'  ',num2str(cpupper(i)),'\n'];
        else
            if cpupper(i) > 0
                textv = [textv,num2str(i),'  Infinity \n'];
            else
                textv = [textv,num2str(i),' -Infinity \n'];
            end
        end
    end
    % i = nums:
    if abs(cpupper(nums)) ~= inf;
        textv = [textv,num2str(nums),'   ',num2str(cpupper(nums)),';\n'];
    else
        if cpupper(nums) > 0
            textv = [textv,num2str(nums),'  Infinity; \n'];
        else
            textv = [textv,num2str(nums),' -Infinity; \n'];
        end
    end
    textv = [textv,'\n'];
    textv = [textv,'param seqlowercp := \n'];
    for i = 1:nums - 1
        if abs(cplower(i)) ~= inf;
            textv = [textv,num2str(i),'  ',num2str(cplower(i)),'\n'];
        else
            if cplower(i) > 0
                textv = [textv,num2str(i),'  Infinity \n'];
            else
                textv = [textv,num2str(i),' -Infinity \n'];
            end
        end
    end
    % i = nums:
    if abs(cplower(nums)) ~= inf;
        textv = [textv,num2str(nums),'   ',num2str(cplower(nums)),';\n'];
    else
        if cplower(nums) > 0
            textv = [textv,num2str(nums),'  Infinity; \n'];
        else
            textv = [textv,num2str(nums),' -Infinity; \n'];
        end
    end
    textv = [textv,'\n'];

end

%%%%%%%%%%%%%%%%%%%%%%%55
% end of cp bounds, beginning of z bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     textv = [textv,'param sequppercp := \n'];
%     for i = 1:nums - 1
%         textv = [textv,num2str(i),'  ',num2str(cpupper(i)),'\n'];
%     end
%     % i = nums:
%     textv = [textv,num2str(nums),'   ',num2str(cpupper(nums)),';\n'];
%     textv = [textv,'\n'];
%      textv = [textv,'param seqlowercp := \n'];
%     for i = 1:nums - 1
%         textv = [textv,num2str(i),'  ',num2str(cplower(i)),'\n'];
%     end
%     % i = nums:
%     textv = [textv,num2str(nums),'   ',num2str(cplower(nums)),';\n'];
%     textv = [textv,'\n'];
    
    % zbounds!!!
    if ((individualopt == 1) || (individualopt == 3)) %individual z bounds
        sizep = size(zlower,1);
        textv = [textv,'# individual bounds for z \n'];
        textv = [textv,'param seqlowerz := \n'];
        for i = 1:sizep - 1
            if (abs(zlower(i)) ~= inf)
                textv = [textv,num2str(i),'  ',num2str(zlower(i)),'\n'];
            else
                if zlower(i) > 0
                    textv = [textv,num2str(i),'  Infinity \n'];
                else
                    textv = [textv,num2str(i),' -Infinity \n'];
                end
            end
        end
        if abs(zlower(sizep)) ~= inf;
            textv = [textv,num2str(sizep),'   ',num2str(zlower(sizep)),';\n'];
        else
            if zlower(sizep) > 0
                textv = [textv,num2str(sizep),'  Infinity; \n'];
            else
                textv = [textv,num2str(sizep),' -Infinity; \n'];
            end
        end
        textv = [textv,'\n'];

        textv = [textv,'param sequpperz := \n'];
        for i = 1:sizep - 1
            if (abs(zupper(i)) ~= inf)
                textv = [textv,num2str(i),'  ',num2str(zupper(i)),'\n'];
            else
                if zupper(i) > 0
                    textv = [textv,num2str(i),'  Infinity \n'];
                else
                    textv = [textv,num2str(i),' -Infinity \n'];
                end
            end
        end
        if abs(zupper(sizep)) ~= inf;
            textv = [textv,num2str(sizep),'   ',num2str(zupper(sizep)),';\n'];
        else
            if zupper(sizep) > 0
                textv = [textv,num2str(sizep),'  Infinity; \n'];
            else
                textv = [textv,num2str(sizep),' -Infinity; \n'];
            end
        end
        textv = [textv,'\n'];

    else %end individualopt == 1

        if zlu == 0
            textv = [textv,'param seqbothz := \n'];
            for i = 1:nums - 1
                if abs(zboth(i)) ~= inf
                    textv = [textv,num2str(i),'  ',num2str(zboth(i)),'\n'];
                else
                    if zboth(i) > 0
                        textv = [textv,num2str(i),'  Infinity \n'];
                    else
                        textv = [textv,num2str(i),' -Infinity \n'];
                    end
                end
            end
            % i = nums:
            if abs(zboth(nums)) ~= inf
                textv = [textv,num2str(nums),'   ',num2str(zboth(nums)),';\n'];
            else
                if zboth(nums) > 0
                    textv = [textv,num2str(nums),'  Infinity; \n'];
                else
                    textv = [textv,num2str(nums),' -Infinity; \n'];
                end
            end

            textv = [textv,'\n'];

        else %% if zlu == 1, we have a zlower and a zupper!

            textv = [textv,'param seqlowerz := \n'];
            for i = 1:nums - 1
                if abs(zlower(i)) ~= inf
                    textv = [textv,num2str(i),'  ',num2str(zlower(i)),'\n'];
                else
                    if zlower(i) > 0
                        textv = [textv,num2str(i),'  Infinity \n'];
                    else
                        textv = [textv,num2str(i),' -Infinity \n'];
                    end
                end
            end
            % i = nums:
            if abs(zlower(nums)) ~= inf
                textv = [textv,num2str(nums),'   ',num2str(zlower(nums)),';\n'];
            else
                if zlower(nums) > 0
                    textv = [textv,num2str(nums),'  Infinity; \n'];
                else
                    textv = [textv,num2str(nums),' -Infinity; \n'];
                end
            end

            textv = [textv,'\n'];

            textv = [textv,'param sequpperz := \n'];
            for i = 1:nums - 1
                if abs(zupper(i)) ~= inf
                    textv = [textv,num2str(i),'  ',num2str(zupper(i)),'\n'];
                else
                    if zupper(i) > 0
                        textv = [textv,num2str(i),'  Infinity \n'];
                    else
                        textv = [textv,num2str(i),' -Infinity \n'];
                    end
                end
            end
            %i = nums:
            if abs(zupper(nums)) ~= inf
                textv = [textv,num2str(nums),'   ',num2str(zupper(nums)),';\n'];
            else
                if zupper(nums) > 0
                    textv = [textv,num2str(nums),'  Infinity; \n'];
                else
                    textv = [textv,num2str(nums),' -Infinity; \n'];
                end
            end

            textv = [textv,'\n'];

        end % end if zlu = 0, else...
    end % end if individual opt = 1 or 3
  
   fprintf(fid,textv);
    
    % let global variable know what type of potparams we currently have
    if zlu == 1     
       potparamswzupperlower = 1;
    else  
        potparamswzupperlower = 0;
    end
    
    
    fclose(fid);
       