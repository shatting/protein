function makepotparamsdotdat(name,pot,nums,numg)
% writes the main data file for the optimization problem
% potential contains ent, R, and mu
% nums is the number of secondary stucture classes (now 27)
% numg is the number of amino acid sequence classes (gamma classes, now 27)
% name is usually potparams.dat
% in ent, mu, R : s's are rows, g's are columns
% clower,cupper,tlower,tupper are the upper and lower bounds for ci, ti,
% dependent on si


global potparamswzupperlower;

fid=fopen(name,'w'); 

    textv = ['# created by makepotparamsdotdat\n'];
    textv=[textv,'param nums:=',num2str(nums),';\n'];
    textv=[textv,'param numg:=',num2str(numg),';\n'];
    textv = [textv,'\n'];
    textv = [textv,'## Entropy\n'];
    textv = [textv,'param ent:'];
    for g=1:numg
        textv = [textv,'   ',num2str(g)];
    end
    textv = [textv,':= \n'];
    for s = 1:nums
        textv = [textv,num2str(s),'   '];
        for g = 1:numg
            x = sysXto10([s g],[nums numg]);
            ent = pot(x).ent;
            if ent ~= inf;
               textv = [textv,num2str(ent),'   '];
            else 
                textv = [textv,'Infinity   '];
            end 
            if g == numg
                if s ~= nums
                    textv = [textv,'\n'];
                else 
                    textv = [textv,';\n'];
                end
            end
        end %end for g
    end %end for s
    
    fprintf(fid,textv);
    
    textv = ['\n'];
    textv = [textv,'## Mean (mu)\n'];
    textv = [textv,'param mu1:'];
    for g=1:numg
        textv = [textv,'   ',num2str(g)];
    end
    textv = [textv,':= \n'];
    for s = 1:nums
        textv = [textv,num2str(s),'   '];
        for g = 1:numg
            x = sysXto10([s g],[nums numg]);
            mu = pot(x).mean;
            mu1 = mu(1);
            if mu1 ~= inf;
               textv = [textv,num2str(mu1),'   '];
            else 
                textv = [textv,'Infinity   '];
            end 
            if g == numg
                if s ~= nums
                    textv = [textv,'\n'];
                else 
                    textv = [textv,';\n'];
                end
            end
        end %end for g
    end %end for s
    textv = [textv,'\n'];
    
    textv = [textv,'param mu2:'];
    for g=1:numg
        textv = [textv,'   ',num2str(g)];
    end
    textv = [textv,':= \n'];
    for s = 1:nums
        textv = [textv,num2str(s),'   '];
        for g = 1:numg
            x = sysXto10([s g],[nums numg]);
            mu = pot(x).mean;
            mu2 = mu(2);  
            if mu2 ~= inf;
               textv = [textv,num2str(mu2),'   '];
            else 
                textv = [textv,'Infinity   '];
            end 
            if g == numg
                if s ~= nums
                    textv = [textv,'\n'];
                else 
                    textv = [textv,';\n'];
                end
            end
        end %end for g
    end %end for s
    textv = [textv,'\n'];
    textv = [textv,'\n'];
    
    textv = [textv,'param mu3:'];
    for g=1:numg
        textv = [textv,'   ',num2str(g)];
    end
    textv = [textv,':= \n'];
    for s = 1:nums
        textv = [textv,num2str(s),'   '];
        for g = 1:numg
            x = sysXto10([s g],[nums numg]);
            mu = pot(x).mean;
            mu3 = mu(3);
            if mu3 ~= inf;
               textv = [textv,num2str(mu3),'   '];
            else 
                textv = [textv,'Infinity   '];
            end 
            if g == numg
                if s ~= nums
                    textv = [textv,'\n'];
                else 
                    textv = [textv,';\n'];
                end
            end
        end %end for g
    end %end for s
    textv = [textv,'\n'];
    
    fprintf(fid,textv);
    textv = ['\n'];
    textv = [textv,'## R \n'];
    textv = [textv,'#at the moment, all of the Rs are triangular - so really do not need all of entries\n'];
     textv = [textv,'param R11:'];
    for g=1:numg
        textv = [textv,'   ',num2str(g)];
    end
    textv = [textv,':= \n'];
    for s = 1:nums
        textv = [textv,num2str(s),'   '];
        for g = 1:numg
            x = sysXto10([s g],[nums numg]);
            r = pot(x).L; %maybe this is L is R^-1, exactly what I think it's supposed to be
            r = r(1,1);
            if r ~= inf;
               textv = [textv,num2str(r),'   '];
            else 
                textv = [textv,'Infinity   '];
            end 
            if g == numg
                if s ~= nums
                    textv = [textv,'\n'];
                else 
                    textv = [textv,';\n'];
                end
            end
        end %end for g
    end %end for s
    textv = [textv,'\n'];
    
%        textv = [textv,'param R12:'];
%     for g=1:numg
%         textv = [textv,'   ',num2str(g)];
%     end
%     textv = [textv,':= \n'];
%     for s = 1:nums
%         textv = [textv,num2str(s),'   '];
%         for g = 1:numg
%             x = sysXto10([s g],[nums numg]);
%             r = pot(x).L; %maybe this is L^-1 or something else -- check later!
%             r = r(1,2);
%             if r ~= inf;
%                textv = [textv,num2str(r),'   '];
%             else 
%                 textv = [textv,'Infinity   '];
%             end 
%             if g == numg
%                 if s ~= nums
%                     textv = [textv,'\n'];
%                 else 
%                     textv = [textv,';\n'];
%                 end
%             end
%         end %end for g
%     end %end for s
%     textv = [textv,'\n'];
%     
%        textv = [textv,'param R13:'];
%     for g=1:numg
%         textv = [textv,'   ',num2str(g)];
%     end
%     textv = [textv,':= \n'];
%     for s = 1:nums
%         textv = [textv,num2str(s),'   '];
%         for g = 1:numg
%             x = sysXto10([s g],[nums numg]);
%             r = pot(x).L; %maybe this is L^-1 or something else -- check later!
%             r = r(1,3);
%             if r ~= inf;
%                textv = [textv,num2str(r),'   '];
%             else 
%                 textv = [textv,'Infinity   '];
%             end 
%             if g == numg
%                 if s ~= nums
%                     textv = [textv,'\n'];
%                 else 
%                     textv = [textv,';\n'];
%                 end
%             end
%         end %end for g
%     end %end for s
    textv = [textv,'\n'];
    
       textv = [textv,'param R21:'];
    for g=1:numg
        textv = [textv,'   ',num2str(g)];
    end
    textv = [textv,':= \n'];
    for s = 1:nums
        textv = [textv,num2str(s),'   '];
        for g = 1:numg
            x = sysXto10([s g],[nums numg]);
            r = pot(x).L; %maybe this is L^-1 or something else -- check later!
            r = r(2,1);
            if r ~= inf;
               textv = [textv,num2str(r),'   '];
            else 
                textv = [textv,'Infinity   '];
            end 
            if g == numg
                if s ~= nums
                    textv = [textv,'\n'];
                else 
                    textv = [textv,';\n'];
                end
            end
        end %end for g
    end %end for s
    textv = [textv,'\n'];
    
       textv = [textv,'param R22:'];
    for g=1:numg
        textv = [textv,'   ',num2str(g)];
    end
    textv = [textv,':= \n'];
    for s = 1:nums
        textv = [textv,num2str(s),'   '];
        for g = 1:numg
            x = sysXto10([s g],[nums numg]);
            r = pot(x).L; %maybe this is L^-1 or something else -- check later!
            r = r(2,2);
            if r ~= inf;
               textv = [textv,num2str(r),'   '];
            else 
                textv = [textv,'Infinity   '];
            end 
            if g == numg
                if s ~= nums
                    textv = [textv,'\n'];
                else 
                    textv = [textv,';\n'];
                end
            end
        end %end for g
    end %end for s
    textv = [textv,'\n'];
    
%     textv = [textv,'param R23:'];
%     for g=1:numg
%         textv = [textv,'   ',num2str(g)];
%     end
%     textv = [textv,':= \n'];
%     for s = 1:nums
%         textv = [textv,num2str(s),'   '];
%         for g = 1:numg
%             x = sysXto10([s g],[nums numg]);
%             r = pot(x).L; %maybe this is L^-1 or something else -- check later!
%             r = r(2,3);
%             if r ~= inf;
%                textv = [textv,num2str(r),'   '];
%             else 
%                 textv = [textv,'Infinity   '];
%             end 
%             if g == numg
%                 if s ~= nums
%                     textv = [textv,'\n'];
%                 else 
%                     textv = [textv,';\n'];
%                 end
%             end
%         end %end for g
%     end %end for s
    textv = [textv,'\n'];
    
       textv = [textv,'param R31:'];
    for g=1:numg
        textv = [textv,'   ',num2str(g)];
    end
    textv = [textv,':= \n'];
    for s = 1:nums
        textv = [textv,num2str(s),'   '];
        for g = 1:numg
            x = sysXto10([s g],[nums numg]);
            r = pot(x).L; %maybe this is L^-1 or something else -- check later!
            r = r(3,1);
            if r ~= inf;
               textv = [textv,num2str(r),'   '];
            else 
                textv = [textv,'Infinity   '];
            end 
            if g == numg
                if s ~= nums
                    textv = [textv,'\n'];
                else 
                    textv = [textv,';\n'];
                end
            end
        end %end for g
    end %end for s
    textv = [textv,'\n'];
    
       textv = [textv,'param R32:'];
    for g=1:numg
        textv = [textv,'   ',num2str(g)];
    end
    textv = [textv,':= \n'];
    for s = 1:nums
        textv = [textv,num2str(s),'   '];
        for g = 1:numg
            x = sysXto10([s g],[nums numg]);
            r = pot(x).L; %maybe this is L^-1 or something else -- check later!
            r = r(3,2);
            if r ~= inf;
               textv = [textv,num2str(r),'   '];
            else 
                textv = [textv,'Infinity   '];
            end 
            if g == numg
                if s ~= nums
                    textv = [textv,'\n'];
                else 
                    textv = [textv,';\n'];
                end
            end
        end %end for g
    end %end for s
    textv = [textv,'\n'];
    
       textv = [textv,'param R33:'];
    for g=1:numg
        textv = [textv,'   ',num2str(g)];
    end
    textv = [textv,':= \n'];
    for s = 1:nums
        textv = [textv,num2str(s),'   '];
        for g = 1:numg
            x = sysXto10([s g],[nums numg]);
            r = pot(x).L; %maybe this is L^-1 or something else -- check later!
            r = r(3,3);
            if r ~= inf;
               textv = [textv,num2str(r),'   '];
            else 
                textv = [textv,'Infinity   '];
            end 
            if g == numg
                if s ~= nums
                    textv = [textv,'\n'];
                else 
                    textv = [textv,';\n'];
                end
            end
        end %end for g
    end %end for s
    textv = [textv,'\n'];
    
    
     
    fprintf(fid,textv);
    
    % let global variable know what type of potparams we currently have
   
    
    fclose(fid);
       
        