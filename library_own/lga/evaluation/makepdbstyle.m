function makepdbstyle(x,seq,fid,name,makestuffup)

sizex = size(x,1); %number of rows of x -- number of amino acids
n = round(sizex/4);
textw=[];
textw=[textw,'MOLECULE ',name,'\n'];
for i = 1:sizex
    aa = seq(i); %will eventually have to change if more atoms are being used
    %model type
    textw=[textw,'ATOM  '];
    fprintf(fid,textw);
    textw=[];
    %atom serial number
    if i < 10
       % textw=[textw,'   ',i];
        fprintf(fid,'    %i',i);
    end
    if ((i > 9) && (i < 100))
       % textw=[textw,'  ',i];
        fprintf(fid,'   %i',i);
    end
    if ((i > 99) && (i<1000))
       % textw=[textw,' ',i];
        fprintf(fid,'  %i',i);
    end
    if (i>1000)
        %textw=[textw,i];
        fprintf(fid,' %i',i);
    end
    %blank space
    textw=[textw,' '];
    %atom name (calpha,c,N,..)
    %for now, just do Ca's
    textw=[textw,' CA '];
    %alternative location index -skip
    textw=[textw,' '];
    %residue name
    if aa < 21
       textw=[textw,aaname(aa)];
    else
        textw=[textw,'   '];
    end
    % blank space plus chain identifier plus residue sequence identifier(make up chains for now)
    if makestuffup == 1
        if i < 10
           textw = [textw,'     ',num2str(i)];
        end
        if ((10 <= i) && (i<100))
            textw = [textw,'    ',num2str(i)];
        end
        if ((100 <= i) && (i<1000))
            textw = [textw,'   ',num2str(i)];
        end
    else
       textw=[textw,'      ']; 
    end
    %code for insertion of residues (skip)
    textw=[textw,' '];
    %blank space
    textw=[textw,'   '];
    
    fprintf(fid,textw);
    textw =[];
    
    %coordinates for x in angstroms
    %for 1-digit before decimal
    if abs(x(i,1)) < 10
        if (x(i,1) >= 0)
           fprintf(fid,'   %.3f',x(i,1));
        else 
           fprintf(fid,'  %.3f',x(i,1));
        end
    end
    %for 2-digits before decimal
    if (( 10 <= abs(x(i,1)) ) && ( abs(x(i,1)) < 100 ))
        if (x(i,1) >= 0)
           fprintf(fid,'  %.3f',x(i,1));
        else 
           fprintf(fid,' %.3f',x(i,1));
        end
    end
    %for 3-digits before decimal
    if (( 100 <= abs(x(i,1)) ) && ( abs(x(i,1)) < 1000 ))
        if (x(i,1) >= 0)
           fprintf(fid,' %.3f',x(i,1));
        else 
           fprintf(fid,'%.3f',x(i,1));
        end
    end
    %for 4digits before decimal
    if ( 1000 <= abs(x(i,1)) ) 
        if (x(i,1) >= 0)
           fprintf(fid,'%.3f',x(i,1));
        else 
           fprintf(fid,'%.2f',x(i,1)); %not enough space for all the decimals
        end
    end
    
    %coordinates for y in angstroms
    %for 1-digit before decimal
    if abs(x(i,2)) < 10
        if (x(i,2) >= 0)
           fprintf(fid,'   %.3f',x(i,2));
        else 
           fprintf(fid,'  %.3f',x(i,2));
        end
    end
    %for 2-digits before decimal
    if (( 10 <= abs(x(i,2)) ) && ( abs(x(i,2)) < 100 ))
        if (x(i,2) >= 0)
           fprintf(fid,'  %.3f',x(i,2));
        else 
           fprintf(fid,' %.3f',x(i,2));
        end
    end
    %for 3-digits before decimal
    if (( 100 <= abs(x(i,2)) ) && ( abs(x(i,2)) < 1000 ))
        if (x(i,2) >= 0)
           fprintf(fid,' %.3f',x(i,2));
        else 
           fprintf(fid,'%.3f',x(i,2));
        end
    end
    %for 4digits before decimal
    if ( 1000 <= abs(x(i,2)) ) 
        if (x(i,2) >= 0)
           fprintf(fid,'%.3f',x(i,2));
        else 
           fprintf(fid,'%.2f',x(i,2)); %not enough space for all the decimals
        end
    end
    
    %coordinates for z in angstroms
    %for 1-digit before decimal
    if abs(x(i,3)) < 10
        if (x(i,3) >= 0)
           fprintf(fid,'   %.3f',x(i,3));
        else 
           fprintf(fid,'  %.3f',x(i,3));
        end
    end
    %for 2-digits before decimal
    if (( 10 <= abs(x(i,3)) ) && ( abs(x(i,3)) < 100 ))
        if (x(i,3) >= 0)
           fprintf(fid,'  %.3f',x(i,3));
        else 
           fprintf(fid,' %.3f',x(i,3));
        end
    end
    %for 3-digits before decimal
    if (( 100 <= abs(x(i,3)) ) && ( abs(x(i,3)) < 1000 ))
        if (x(i,3) >= 0)
           fprintf(fid,' %.3f',x(i,3));
        else 
           fprintf(fid,'%.3f',x(i,3));
        end
    end
    %for 4digits before decimal
    if ( 1000 <= abs(x(i,3)) ) 
        if (x(i,3) >= 0)
           fprintf(fid,'%.3f',x(i,3));
        else 
           fprintf(fid,'%.2f',x(i,3)); %not enough space for all the decimals
        end
    end
    
    %occupancy
    if makestuffup == 1
        textw=['  1.00']; %make something up for maxsprout
    else 
        textw = ['      '];
    end
    fprintf(fid,textw);
    textw=[];
    %temperature factor
    if makestuffup == 1
        textw=[textw,'  0.00'];
    else
        fprintf(fid,'      ');
    end
      %textw=[textw,'      '];
%     if makestuffup == 1
%        %blank space for maxsprout
%        textw=[textw,'      '];
%        Crn whatever for maxsprout
%        textw=[textw,'1CRN'];
%        line number for maxsprout
%        if i + 69 < 100
%           textw=[textw,'  ',num2str(i+69)];
%        end
%        if ((i + 69 >= 100) && (i + 69)<1000)
%           textw=[textw,' ',num2str(i+69)];
%        end
%        if i + 69 > 1000
%           textw=[textw,num2str(i+69)];
%        end
%    else
%        %blank space - for not maxsprout
%        textw=[textw,'          '];
%        %element symbol, right-justified
%         if aa < 21
%           [name,letter]=aaname(aa);
%           textw=[textw,' ',letter];
%         else
%             textw=[textw,'  '];
%         end
%     end
%     %charge on the atom
     textw=[textw,'  \n'];
     fprintf(fid,textw);
     textw=[];
%     
end %end i=1:sizex
textw=[textw,'END\n'];
fprintf(fid,textw);
textw=[];
