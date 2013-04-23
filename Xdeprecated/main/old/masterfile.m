global proteinroot;

OPTIONS_CLASS.max_iterations = 7;
OPTIONS_CLASS.plots = 0;    % do plots?
OPTIONS_CLASS.manual_t = 0; % use manual initial t-classes? (memory or file)
OPTIONS_CLASS.ask = 0;      % ask each iteration of classification?

if ~exist('data','var')
    load RT071127;
end

if ~exist('pot') % || exist('extlookup') == 0,
    mf = [proteinroot,filesep,'masterdata.mat'];
    if exist(mf,'file') && ask('load from masterfile?'),
        disp('loading masterdata from file.');
        load(mf);
        disp('done.');
    else
        %fourfrag_suff_linear,pot, cl(:,i), confus, freq
        
        if ~exist('pot'),
            getpotential;
        end
        
        
        if ~exist('extlookup') && 0, %no extlookup as of now(24.3)
            getcatfeatures;
        end
        
        lookup = createlookup(fourfrag_suff_linear, cl);

        vsort = vgamma(lookup(:,:,:,:,end),pot,data);
        
        %[extlookup,extpot,confusion]

        %TODO: rebuild potentials from lookup
        %not done because we use lookup for now just for unknown fragments
        
        save(mf, 'lookup', 'pot', 'cl','vsort');  
        
    end
end

p = 137; %try 137,323,194,191,190,188,175,169 -- all smallish

createampl('tpltest', data, p, lookup(:,:,:,:,end), pot, 0.99, 0.05, 0,1,1,0, vsort,1,0);
%    (probname, data, p, lookup, pot, alpha, beta, initials, avginitials, newconstraint,makefeasible, Vsort,algnum,maxit )

if (ispc)
    !optimize
else
    !./optimize.sh
end

getandcompare(p,lookup,pot,data);

% fragclasses = classseq( data{p}.seq, lookup );
% 
% fragclasses = fragclasses(:,end:-1:1); % reverse iteration order: final classes first
% 
% figure;
% maximize_fig;
% title 'predicted fold';
% ribbon3dnice(bond2coords(geometry2bondn(fromampl('tpltest'))), fragclasses, data{p});
% 
% figure;
% maximize_fig;
% title 'actual fold';
% ribbon3dnice(bond2coords(data{345}.bond), fragclasses, data{p});