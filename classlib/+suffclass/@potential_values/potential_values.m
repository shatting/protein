classdef potential_values < handle
%  suffclass.potential_values < handle
%   Vectorized class prediction and potential evaluation. 
%
% PUBLIC PROPERTIES
%   .potential  suffclass.potential handle
%   .options    
%       .ent_weight (default 1)
%       .handleequalpotentials
%           false, "warn", "error", "unclassified" or "randomize" 
%           (TODO: randomize partially implemented)
%       .use_potential_entropy (1)
%   .priorsORentropies      saved value from last .calcvals invocation
%   .ispriors               saved value from last .calcvals invocation
%
% PUBLIC TRANSIENT PROPERTIES
%   .Vdiff
%   .cl(l,1)        predicted class for lth data vetor
%   .vmin(l,1)      minimum potential value
%   .V(l,cl)        potential value of lth vector under potential cl
%   .prob(l,cl)     probability that lth vector is in class cl
%
% CONSTRUCTORS
%   this = potential_values(@potential,options)
%
% PUBLIC METHODS    
%   cl = .calcvals(data,priorsORentropies,ispriors) (was suff_data2v)
%       Calculate class assignment, potential values and conditional
%       probabilities for data under potential pot. This is the vectorized
%       version of suff_data2cl.
%       
%       INPUT
%       data(l,:)               l-th data vector
%       priorsORentropies(cl,1) (optional) input prior or entropy
%                               (log priors) of class cl. if not given, 
%                               then uniform assumed
%       ispriors

    
    properties
        options
        potential
        
        priorsORentropies
        ispriors
    end
    
    properties (Transient)                
        cl
        vmin
        V
        prob
        Vdiff
    end
    
    methods        
        function this = potential_values(potential,options)
            this.potential = potential;
            this.options = suffclass.potential_values.defaultoptions(options);
        end
                     
        function cl = calcvals(this,data,priorsORentropies,ispriors)
            
            suff = this.potential.suff;
            nnum = suff.nnum;
            ncat = suff.ncat;
                        
            
            
            if nargin>=3 && ~isempty(priorsORentropies),
                if length(priorsORentropies)~=suff.ncl,
                    suff.ncl,size(priorsORentropies)
                    error('prior length must equal number of classes');
                end;
                if (nargin==3), error('calcvals signature changed'); end
            end;
            
            for i=1:ncat,
                ii=suff.cat(i);
                if min(data(:,ii))<=0,
                    error('categorical variables must be positive')
                else
                    [dmax,ind]=max(data(:,ii));
                    if dmax>abs(suff.typ(ii)),
                        ind,ii,dmax
                        error('data(ind,ii) exceeds number of levels')
                    end;
                end;
            end;
            
            N=size(data,1);
            on=ones(N,1);
            numpot = this.potential.numpot;
            catpot = this.potential.catpot;
            ncl = suff.ncl;
            
            for cli=ncl:-1:1,
                % decide 
                if nargin==2 % use uniform prior
                    entcl=zeros(N,1);
                elseif sum(priorsORentropies) ~= 0 && ~this.options.use_potential_entropy,                    
                    if ~ispriors, % use input entropies
                        entcl=zeros(N,1) + this.options.ent_weight*priorsORentropies(cli);
                    else % use input probabilities
                        entcl=zeros(N,1) + this.options.ent_weight*log(priorsORentropies(cli)+realmin);
                    end
                else % use entropies from potentials
                    %TODO: either consensus classification or figure out a way to add
                    %cat and num potentials without adding prior 2 times                    
                    catpotent_i = 0;
                    numpotent_i = 0;
                    if (~isempty(catpot)), catpotent_i = catpot.ent(cli); end
                    if (~isempty(numpot)), numpotent_i = 0.5*log(det(numpot(cli).cov)); end; %numpot(cli).ent; end
                    
                    entcl = zeros(N,1) + this.options.ent_weight*(catpotent_i + numpotent_i);
                end
                
                numpotcl = zeros(N,1);
                if nnum>0,
                    row=numpot(cli).mean';
                    diff=data(:,suff.num)-row(on,:);
                    diffR=diff*numpot(cli).L';
                    for k=1:nnum,
                        numpotcl=numpotcl+diffR(:,k).*diffR(:,k);
                    end;
                    numpotcl = numpotcl/2 +  length(suff.num)/2*log(2*pi);
                end;
                
                catpotcl = zeros(1,N);
                if ncat>0,                    
                    for i=1:ncat,
                        catpotcl=catpotcl+catpot.pot(i,data(:,suff.cat(i)),cli);
                    end                    
                end
                potcl = entcl + catpotcl' + numpotcl;
                V(:,cli)=potcl;
                
            end
            [this.vmin,cl]=min(V,[],2); %TODO: equal potentials: randomize class assignment?
            this.cl = cl;
            this.V = V;
            
            if this.options.handleequalpotentials
                Vdiff = sort(abs(V - repmat(this.vmin,1,ncl)),2); % sort potential differences; then the first column will always be zero
                Vdiff = min(Vdiff(:,2:end),[],2);
                Vequal = Vdiff < 1e-10;
                if any(Vequal),
                    warning(sprintf('%i (%d%%) of data points had two or more classes with the same potential value.',sum(Vequal),sum(Vequal)/length(Vequal)*100));
                end
                this.Vdiff = Vdiff;
            end
            
            %if nargout<4, return; end;
            % get probability table
            for cli=ncl:-1:1,
                epot(:,cli)=exp(this.vmin-V(:,cli));
            end;
            epotsum = sum(epot,2);
            this.prob=epot./epotsum(:,ones(ncl,1));
            
        end
    end
    
    methods (Static)
        function options = defaultoptions(options)
            
            if (nargin==0)
                options = struct;
            end
            
            if ~isfield(options,'ent_weight'), options.ent_weight = 1; end
            if ~isfield(options,'use_potential_entropy'), options.use_potential_entropy = 0; end
            if ~isfield(options,'handleequalpotentials'), options.handleequalpotentials= 0; end
            
        end
    end
    
end

