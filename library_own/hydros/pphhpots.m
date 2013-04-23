function pot = pphhpots(feats)
% uses a covariance model to make a potential for the hydrophilic
% and hydrophobic aa's. based on feat14 (centeredness feature) and 
% feat15 (edginess feature), which are the 14th and 15th rows in Y, 
% and are stored separately in hydrofeats.mat
% not based on feat15 anymore, the formula is too complicated for me to 
% type into ampl form now. I'll do it later if it seems useful
% so feats is now based on just feat14
% I divided ent by 1000 bc it was so big that the other numbers were losing
% their impact

% so load hydrofeats.mat first to get feats

% how hydrofeats was made (so I know how to quickly change this to using
% more than one feature)
% load allfeats
% feats = cell(20,1);
% for a = 1:20
%     feats{a} = Y{a}(14:15,:);
%     clear Y{a};
% end
% clear Y;


pot=struct('ent',0,'L',0,'mu',0);

for a = 1:20
    centerednessf = feats{a}(:,14);
    % first need to center the feature matrix
    n = size(feats{a},2);
    S = sum(feats{a},2); % S is a (2*1) vector
    ybar = S/n;
    %an = ones(n);
    ybarm = repmat(ybar,1,n); % [ybar,ybar,ybar,...]
    F = feats{a} - ybarm; % the centered matrix!
    clear feats{a};
    
    %now the covariance model
    [Q,R] = qr(F',0); %Q is n*2, R is 2*2
    enta = log(det(R)^2)/10000;
    
    pot(a).ent = enta; % divide by 10000 bc its so much bigger than the other numbers, that they lose their significance
    pot(a).L = (R^(-1))';
    pot(a).mu = ybar;
    
end % for all a

Readme=[
'made by pphhpots.m, based on feat14 for each a           '
'pot(a).ent = gamma = constant                            '
'pot(a).L = 2*2 lower triangular matrix                   '
'pot(a).mu = ybar = 2*1 vector                            '
'divided ent by 10000 bc it was outpowering the others    '
];

save hydropots pot Readme;

if 1
% make a polarities matrix to compare pot values with hydrophob, hydrophil
polarities = zeros(1,20);
for i = 1:20
[name,letter,polarity,all]=aaname(i);
if polarity == 'PP'
polarities(i) = 1;
elseif polarity == 'PO'
polarities(i) = 1;
elseif polarity == 'HH'
polarities(i) = 2;
elseif polarity == 'Ho'
polarities(i) = 2;
elseif polarity == 'HO'
polarities(i) = 2;
elseif polarity == 'P1'
polarities(i) = 1;
else
polarities(i) = 0;
end
end

if (size(pot(1).L,1) == 1 && (size(pot(1).L,2) == 1))
 maybesuccess = zeros(20,4);
   for i = 1:20
    maybesuccess(i,1) = polarities(i);
    maybesuccess(i,2) = pot(i).mu;
    maybesuccess(i,3) = pot(i).L;
    maybesuccess(i,4) = pot(i).ent;
   end
else
  maybesuccess = zeros(20,4);
   for i = 1:20
    maybesuccess(i,1) = polarities(i);
    maybesuccess(i,2) = pot(i).mu(1);
    maybesuccess(i,3) = pot(i).mu(2);
    maybesuccess(i,4) = pot(i).ent;
   end
end
disp('pol (1: PP  2: HH)            mu1             mu2 or L (depending on size)      ent');
  maybesuccess

end % end if 0


    
    