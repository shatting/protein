function lgacoords = readoutLGA(filename)
% reads coordinates of superpositioned template file from LGA
% careful: LGA rotates the template, not the model!!!!
% new: checks for error: 

lga = fopen(filename);

%first find number of aa's
% given in line: REMARK   SUMMARY:   33
isline = 0;

while isline == 0
    thisline = fgetl(lga); %reads in a line of the file
    if size(thisline,2) < 3
        continue;
    end
    linecontains = thisline(10:11);  %reads 11th and 12th letters in this line
    if min(linecontains == 'SU') == 1
        isline = 1;
    end
    if size(thisline,2) >= 8
        if min(thisline(1:8) == '# ERROR!') == 1
            lgacoords = 0;
            return
        end
    end
end

% get protein size
linecontains = str2num(thisline(18:end));
naa = linecontains(1);

%now get calpha prediction
lgacoords = zeros(naa,3);

% find beginning of ATOM section, where coords are listed
isline = 0;
while isline == 0
    thisline = fgetl(lga); %reads in a line of the file
    if size(thisline,2) < 3
        continue;
    end
    linecontains = thisline(1:4);  %reads 11th and 12th letters in this line
    if min(linecontains == 'ATOM') == 1
        isline = 1;
    end
    if size(thisline,2) >= 8
        if min(thisline(1:8) == '# ERROR!') == 1
            lgacoords = 0;
            return
        end
    end
end

% here come the coords:
linecontains = str2num(thisline(28:end));
lgacoords(1,:) = linecontains(1:3);
for i = 2:naa
    thisline = fgetl(lga);
    linecontains = str2num(thisline(28:end));
    lgacoords(i,:) = linecontains(1:3);
end

fclose(lga);