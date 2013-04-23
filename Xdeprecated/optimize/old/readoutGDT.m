function [GDT,rmsd] = readoutGDT(filename)
% reads the GDT_TS value out of the LGA output
% new: checks for ERROR: # ERROR!

lga = fopen(filename);

isLGA = 0;


%find LGA block, since GDT_TS score is directly below it
while (isLGA == 0)
    thisline = fgetl(lga); %reads in a line of the file
    if size(thisline,2) < 3
        continue;
    end
    linecontains = thisline(1:3);  %reads first 3 letters in this line
    if min(linecontains == 'LGA') == 1
        isLGA = 1;
    end
    if size(thisline,2) >= 8
        if min(thisline(1:8) == '# ERROR!') == 1
            GDT = 0;
            rmsd = inf;
            return
        end
    end
end

% now find end of LGA block
while isLGA == 1
   thisline = fgetl(lga); %reads in a line of the file
   if size(thisline,2) < 3
        continue;
    end
    linecontains = thisline(1);  %reads first 3 letters in this file
    if min(linecontains ~= 'L') == 1
        isLGA = 0;
    end
    if size(thisline,2) >= 8
        if min(thisline(1:8) == '# ERROR!') == 1
            GDT = 0;
            rmsd = inf;
            return
        end
    end
end

%now skip two lines to get to GDT line
%thisline = fgetl(lga);
%thisline = fgetl(lga);

%now I'm at the right line!
thisline = fgetl(lga); %reads in a line of the file
linecontains = str2num(thisline(14:end));  %skip first 14 spaces, bc they say 'SUMMARY(GDT)', and I want an array of numbers
GDT = linecontains(6);  %GDT_TS should be the sixth item
rmsd = linecontains(5);

fclose(lga);
%close(filename);