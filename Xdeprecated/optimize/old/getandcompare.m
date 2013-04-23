function [GDT,LGArmse,rmse] = getandcompare(p,data,probname,newclasses,confus,class_info,s)
% if bondls == 1, fromampl gets information about all the variables used in
% ampl to determine the bond lengths

%if nargin < 5
%    bondls = 0;
%end
if nargin < 7
  if newclasses == 0
     s = seq2sclasses(data{p}.bond);
  else
     s = new_deal_hefg_interface(class_info,data{p}.bond,data{p}.seq);
  end
end
%gamma = seq2gammaclasses(data{p}.seq,class_info);
%fragclasses = fragclasses(:,end:-1:1); % reverse iteration order: final classes first

%get folds:
prez = fromampl(probname,s,newclasses,1); %gives prez.zo, .co, .alphaso
if newclasses == 0
   pregeom = geometryz2geometry(prez);
   pre = geometry2bondn(pregeom);%*3.8;
else
   prez.tqu = prez.zo;
   prez.tcl = s;
   prez.cpo = prez.co(end);
   prez.co = prez.co(1:end-1); 
   pre = geo2bond(prez);
end
% fromampl now expects just a c and a z, no cp

%old = double(data{p}.bond);
%old = rot(old); %rotation is new!!
%old = old./repmat(sqrt(sum(old.^2,2)),1,3);
old = rotatematrix(data{p}.bond);

disp('twist with Q and shift by q to make prediction from original configuration:')
[ Q, q, rmse ] = register(old,pre)
%finds optimal superposition of pre onto old

regold = old*Q; %not quite sure what to do with q...
% for i = 1:size(regold,1)
%     regold(i,:) = regold(i,:) + q;
% end
%compare3d(old,pre)
%compare3d(regold,pre)
%rot here is new!!! (took it away: it was calphas(rot(double(data{p}.bond))
%3.8* here is new!
%makepdbform(1,3.8*calphas(pre),mirrorx(3.8*calphas(old)),data{p}.seq);

[GDT,LGArmse,lgacoords] = dolga(pre,data{p});

% cals = pre; %calphas(pre);
% mat1 = bond2mat(data{p});
% mat2 = mat1;
% for i = 1:size(cals,1)
%     mat2.atomr(i).x = cals(i,1);
%     mat2.atomr(i).y = cals(i,2);
%     mat2.atomr(i).z = cals(i,3);
% end
% 
% % get into the correct directory
% % I have no idea how to do this properly
% % see if I am at home:
% global proteinroot;
% athome = 0;
% if min(proteinroot == 'C:\Users\Owner\prot')
%     athome = 1;
% end
% if athome == 1
%     cd C:\Users\Owner\prot\ginny
% end
% 
% filename = writelgainput(mat1,'template',mat2,'model','firsttry');
% delete 'lgaoutput';
% diary('lgaoutput');
% runlga(filename,'-3 -o2'); %options are supposed to come in the ''. want '-3' (and maybe -atom:CA -o1', but I think they're the defaults)
% diary off
% % get GDT
% filename = 'lgaoutput';
% [GDT,LGArmsd] = readoutGDT(filename)
% % gdt LGA's prediction of best superposition of template onto model
% delete 'firsttry'
% cd ..
% movefile('stephan\lga\lga_bin\firsttry.pdb','ginny\firsttry')
% cd ginny
% filename = 'firsttry';
% lgacoords = readoutLGA(filename) 
% fclose('all');

if min(size(lgacoords) == [1 1]) == 1
    disp('LGA found error');
    compare3d(regold,calphas(pre))
else
    compare3d(lgacoords,calphas(pre),1) % 1 at end indicates that already in calpha coordinates
end


% 
% figure;
% %maximize_fig;
% title 'actual fold';
% ribbon3dnice(bond2coords(regold), fragclasses, data{p});
display('first fig is target fold');
% 
% figure;
% %maximize_fig;
% title 'predicted fold';
% ribbon3dnice(bond2coords(pre), fragclasses, data{p});
display('second fig is predicted fold');


%disp('potential value at original configuration:') %this isn't important,
%and needs to be updated 
%V = findV(p,data,confus,class_info)

% disp('Vhydro value after optimization:')
% prez.Vhydro
% display('Vhydro value of original configuration:')
% Vhydro = hydropotential(data,p)

%V = proteinpotential(data{p},lookup(:,:,:,:,end),pot)

%disp('to load bond vectors and calphas coordinates, try');
%disp('fromampl(tpltest,1) (put tpltest in quotations)');
%disp('bondlengthconstraintsfeasibilitytester');

%disp('new thing to do: switch ccs (in createampl, tpltest.mod\n');
%disp('to half squared distances\n');
