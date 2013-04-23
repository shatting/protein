function createproteinfile(s,gamma,probname,p,actualcc,hp,newclasses,thfrags,n3f,thrmax,thrfr,Vfrag,pair2cl)
% ex. createproteinfile(s,gamma,[ampldir,'newprob','_protein.dat'],323);
% creates the protein.dat file for the protein p
% contains length of sequences (n), the s sequence, and the gamma sequence for the protein p
% hp params are in here for now, since I don't always use them

if nargin < 5
    actualcc = 0;
end
if nargin < 6
    hp = 0;
end
if (hp > 0)
    load RT071127.mat
end
if nargin < 7
    thfrags = 0;
end
if nargin < 12
    pair2cl = 0;
end
%if thfrags > 0
%    load thfragsstuff; %loads file tf.pair2cl, tf.potc
%end

n = length(s); % = length of protein - 3 (number of 4-fragments) = length of s seq, gamma seq

fid=fopen(probname,'w');

textv = ['#created by createproteinfile\n'];

textv = [textv,'#protein: ',num2str(p),';\n'];
textv = [textv,'# n is the length of the gamma and s sequences\n'];
textv = [textv,'#param n:=',num2str(n),';\n'];
textv = [textv,'# cis is the number of ci"s in the protein geometry\n'];
textv = [textv,'#param cis:=',num2str(n+1),';\n'];
%cc is inverse sum of close contacts of actual protein (half of the
%contacts, bc symmetric)
textv = [textv,'## cc is used for the close contact constraints\n'];
if actualcc > 0
  textv = [textv,'param cc:=',num2str(actualcc),';\n'];
  if hp > 0
     pol=[1.1417e+003,-2.8793e+003, 1.1208e+004];
     fit= @(x) pol(1).*x.^(2/3)+pol(2).*x.^(1/3)+pol(3);  
     textv = [textv,'param fitp:=',num2str(fit(n+2)),';\n'];
     %the primary sequence
     textv = [textv,'#the primary sequence\n'];
     aaseq = data{p}.seq;
     %clear data;
     textv = [textv,'param aaseq:=\n'];
     for j=1:(n+2),
        textv = [textv,num2str(j),'   ',num2str(aaseq(j)),'\n'];
     end
     textv = [textv,num2str(n+3),'   ',num2str(aaseq(n+3)),';'];
     textv = [textv,'\n'];
     % upper and lower bounds for Vhy (come from hydrobylength, very basic (max and min Vhy based on protein length) 
     minh = @(x) 2.2763*x -101.3816;
     maxh = @(x)  2.0601*x + 89.0545;
     naa = size(aaseq,1);
     textv = [textv,'param hydrodown:= ',num2str(minh(naa)),'; \n'];
     textv = [textv,'param hydroup:= ',num2str(maxh(naa)),'; \n'];
     
  end %end if hp > 0
end % end if actualcc > 0
textv = [textv,'\n'];

fprintf(fid,textv);
textv = [];
textv = [textv,'#the secondary structure class sequence \n'];
if newclasses == 0
    textv = [textv,'param seq:=\n'];
    for j=1:(n-1),
        textv = [textv,num2str(j),'   ',num2str(s(j)),'\n'];
    end
    textv = [textv,num2str(n),'   ',num2str(s(n)),';'];
    textv = [textv,'\n'];
end


textv = [textv,'\n'];
textv = [textv,'#the primary structure class sequence \n'];
textv = [textv,'param gseq:=\n'];
for j=1:(n-1)
    textv = [textv,num2str(j),'   ',num2str(gamma(j)),'\n'];
end
textv = [textv,num2str(n),'   ',num2str(gamma(n)),';'];

textv = [textv,'\n'];

global alphas;
if actualcc > 0
    if newclasses == 0
    alpha = [alphas(:,1:2),alphas(:,3:5)*pi/100];
    cz = [alpha(:,1:2),cos(alpha(:,3:5))];
    sz = [alpha(:,1:2),sin(alpha(:,3:5))];
    y = sys10toX(s,[3 3 3]);
    si = zeros(size(s));
    ci = zeros(size(s));
   % alphao = zeros(size(s));
    for i = 1:size(s,2)
       si(i)  = sz(y(i,2) + (y(i,1) - 1)*3,y(i,3) + 2);
       ci(i) =  cz(y(i,2) + (y(i,1) - 1)*3,y(i,3) + 2);
    %   alphao(i) = alpha(y(i,2) + (y(i,1) - 1)*3, y(i,3) + 2);
    end
    
    textv = [textv,'\n'];
%     textv = [textv,'param alph:=\n'];
%     for j=1:(n-1)
%         textv = [textv,num2str(j),'   ',num2str(alphao(j)),'\n'];
%     end
%     textv = [textv,num2str(n),'   ',num2str(alphao(n)),';'];
%     textv = [textv,'\n'];
%     textv = [textv,'\n'];
   textv = [textv,'## cz and sz are cos, sin(alpha_i) \n'];
   textv = [textv,'## alpha_i is the midpoint of the range for phi_i \n'];
   textv = [textv,'## range for phi_i is dependent on s class  \n'];
   textv = [textv,'## so cz, sz are independent of the value of z, but \n'];
   textv = [textv,'## dependent on the s class - perfect! \n'];
%     textv = [textv,'param cz:=\n'];
%     for j=1:(n-1)
%         textv = [textv,num2str(j),'   ',num2str(ci(j)),'\n'];
%     end
%     textv = [textv,num2str(n),'   ',num2str(ci(n)),';'];
%     textv = [textv,'\n'];
% 
%     textv = [textv,'\n'];
%     textv = [textv,'param sz:=\n'];
%     for j=1:(n-1)
%         textv = [textv,num2str(j),'   ',num2str(si(j)),'\n'];
%     end
%     textv = [textv,num2str(n),'   ',num2str(si(n)),';'];
%     textv = [textv,'\n'];
    end   % ends if newclasses == 0

    
    if hp > 0
        load hydropots;
        %mu
        textv = [textv,'\n'];
        textv = [textv,'param hpmu:=\n'];
        for a=1:19
           textv = [textv,num2str(a),'   ',num2str(pot(a).mu),'\n'];
        end
        textv = [textv,num2str(20),'   ',num2str(pot(20).mu),';'];
        textv = [textv,'\n'];
        textv = [textv,'\n'];
        textv = [textv,'param hpL:=\n'];
        for a=1:19
           textv = [textv,num2str(a),'   ',num2str(pot(a).L),'\n'];
        end
        textv = [textv,num2str(20),'   ',num2str(pot(20).L),';'];
        textv = [textv,'\n'];
        textv = [textv,'\n'];
        textv = [textv,'param hpent:=\n'];
        for a=1:19
           textv = [textv,num2str(a),'   ',num2str(pot(a).ent),'\n'];
        end
        textv = [textv,num2str(20),'   ',num2str(pot(20).ent),';'];
        textv = [textv,'\n'];
        clear pot;
    end % end if hp > 0
    
end % end if actualcc > 0

%three frags parameters! (pot params are in ft.dat, these are
%parameters specific to the protein

if thfrags > 0
%    [threefrags,Vfrag] = get3fragpot(p,data,tf.pair2cl,tf.potc); 
    n3fcl = thrmax; %max(max(tf.pair2cl));
    threefrags = thrfr;
    %n3f = size(threefrags,1);
   
    textv = [textv,'   \n'];
    textv = [textv,'# 3-frags parameters specific to this protein \n'];
    textv = [textv,'param n3fclass := ',num2str(n3fcl),';  # not actually specific to this protein \n'];
    textv = [textv,'#param n3f := ',num2str(n3f),';  #number of 3-fragments in this protein \n'];
    textv = [textv,'   \n'];
    
   textv = [textv,'param thfrags: 1 2 :=   \n'];
    for i = 1:n3f-1
        textv = [textv,num2str(i),'  ',num2str(threefrags(i,1)),'   ',num2str(threefrags(i,2)),' \n'];
    end
   textv = [textv,num2str(n3f),'  ',num2str(threefrags(n3f,1)),'   ',num2str(threefrags(n3f,2)),'; \n'];
    textv = [textv,'   \n'];
    
    textv = [textv,'   \n'];
    textv = [textv,'param thfragsclass:=  \n'];
    for i = 1:n3f - 1
        textv = [textv,num2str(i),'  ',num2str(threefrags(i,5)),' \n'];
    end
    textv = [textv,num2str(n3f),'  ',num2str(threefrags(n3f,5)),'; \n'];
    textv = [textv,'   \n'];
    
%     textv = [textv,'param pair2cl:  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  20:= \n   '];
%     for i = 1:20
%         textv = [textv,num2str(i),'   '];
%         for j = 1:19 
%             textv = [textv,num2str(pair2cl(i,j)),'   '];
%         end
%         if i < 20
%             textv = [textv,num2str(pair2cl(i,20)),'   \n'];
%         else
%             textv = [textv,num2str(pair2cl(i,20)),';   \n'];
%         end
%     end
          
    textv = [textv,'   \n'];        
    textv = [textv,'#param rmaxsq = # have to figure out what the input of the weight function is   \n'];
    textv = [textv,'   \n'];
    textv = [textv,'param thfragsmax = ',num2str(Vfrag),'; #the actual value of thfragspot for this protein   \n'];
    textv = [textv,'   \n'];
    textv = [textv,'   \n'];
    textv = [textv,'   \n'];
    textv = [textv,'   \n'];

end % end if thfrags
    
fprintf(fid,textv);
fclose(fid);

%edit(filename);