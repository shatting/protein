function [ data ] = letter_recognition()

nattrs = 17;

data = zeros(1,nattrs);

fid=fopen('letter-recognition.data');
l = 1;
while 1
    tline = fgetl(fid);    
    if ~ischar(tline) || strcmp(strtok(tline,','),'') ,   break,   end
    
    cll = regexp(tline, ',', 'split');    
    data(l,end) = int32(cll{1}) - int32('A') + 1;
    data(l,1:end-1) = str2num(tline(3:end));
    %disp(tline);
    %disp(data(l,:));
    l = l+1;
end
fclose(fid);

end
