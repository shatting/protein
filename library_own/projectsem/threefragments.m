if exist('data')~=1, 
  tic
  clear
  load RT.mat 
  disp('data loaded')
  showtime
end;

tic

nbin=15;
histogram = zeros(20,20,20,nbin);

% %a=sum(histogram,4);
% %mean(a(:)); = 42.0380
% basecellsize = 40;
% angles = cell(20,20,20,basecellsize);
% anglesinfo = cell(20,20,20);

numproteins = length(data);
nonemptycount=0;

for pr=1:numproteins,
    
    if (mod(pr,500) == 0)
        disp(pr);
    end
    
    % arbeitsvariablen
    prcell = data{pr};
    seq = prcell.seq;
    numaa = length(seq);    
    bond = double(prcell.bond);

    % bonds normalisieren
    lengths = sqrt(sum(bond.*bond,2));
    bondn = bond./[lengths lengths lengths];
    
    % gramn = bondn*bondn';
    
    % datenvektoren
    c = zeros(numaa-1,1); % cos(alpha_i)
    s = zeros(numaa-1,1); % sin(alpha_i)
    gi2 = zeros(numaa-3,1); % tilde{G}_{i,i+2}
    t = zeros(numaa-2,1); % t=cos(beta_i)
    tipr = zeros(numaa-2,1); % t'=sin(beta_i)
    
    for i=1:numaa-2,
        c(i) = bondn(i,:)*bondn(i+1,:)';
        %if (i<numaa-2)
            %gi2(i) = bondn(i,:)*bondn(i+2,:)';
        %end
    end
        
    %s = sqrt(1-c.^2);
    
    %for i=1:numaa-3,
       %t(i) = (c(i)*c(i+1)-gi2(i))/(s(i)*s(i+1));
       %tipr(i) = det([bondn(i,:);bondn(i+1,:);bondn(i+2,:)])/(s(i)*s(i+1));
    %end

    
    % AUSWAHL
    daten = c;       
    datenmin=-1;
    datenmax=1;
    
    % bin setup
    binlen=(datenmax-datenmin)/nbin;    

    for i=1:numaa-2,
        aa1 = seq(i);
        aa2 = seq(i+1);
        aa3 = seq(i+2);
        if (aa1 > 20 || aa2>20 || aa3 > 20)
            continue;
        end                        
                       
        bin=round((daten(i)-datenmin)/binlen+0.5);
        
        histogram(aa1,aa2,aa3,bin) = histogram(aa1,aa2,aa3,bin)+1 ;
        
%         if ~isfield(anglesinfo{aa1,aa2,aa3},'count')
%             anglesinfo{aa1,aa2,aa3}.count = 0;
%             anglesinfo{aa1,aa2,aa3}.cellsize = basecellsize;
%         end
%         
%         count = anglesinfo{aa1,aa2,aa3}.count;
%         cellsize = anglesinfo{aa1,aa2,aa3}.cellsize;
%         
%         if (count == cellsize )
%             angles{aa1,aa2,aa3} = [angles{aa1,aa2,aa3};cell(cellsize,1)];
%             anglesinfo{aa1,aa2,aa3}.cellsize = 2*cellsize;
%         end
%         
%         angles{aa1,aa2,aa3,count + 1} = gramn(i,i+1);
%         anglesinfo{aa1,aa2,aa3}.count = count + 1;
        
    end
    
end

t = toc
showtime
disp(sprintf('proteine/s: %f',numproteins/t));
