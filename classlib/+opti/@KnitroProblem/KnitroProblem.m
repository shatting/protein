classdef KnitroProblem < handle
    % opti.KnitroProblem < handle
    %   Encapsulates data and methods for one knitro run.
    
    properties
        sgseq
        %optidata                
        fevaldata
        avglen = 38; %TODO: more elaborate
        savexhistory
        options
        sclassifier
        
        havestarted
        haverun
        
        aaseq
        sseq
        gseq
        nfrags
        dimx
        ind_c
        ind_z
        fobj
        LB
        UB
        X0
        A
        B
        Aeq
        Beq
        NONLINCON
        
        geompred
        FVAL
        EXITFLAG
        OUTPUT
        LAMBDA
        
        xhistory
        xhistcoords
        predcoords
        preddataset
        knitrotime
    end
    
    methods
        
        function obj = KnitroProblem(aaseq, sseq, gseq, optidata, options)
            if (nargin < 4)
                options = struct;                
            end
            % set default options no, we want errors if no options are set
            %if ~isfield(options,'savexhistory'), options.savexhistory = 0; end
            %if ~isfield(options,'X0'), options.X0 = []; end
            %if ~isfield(options,'useX0'), options.useX0 = 0; end            
            
            obj.options = options;    
            obj.sclassifier = optidata.sclass.classifier;
            obj.sseq = sseq;
            obj.gseq = gseq;
            obj.aaseq = aaseq;
            obj.sgseq = optidata.sgclass.classifier.getcombicl(sseq, gseq);
            
            obj.preparebounds(optidata);            
        end
        
        function savexcallback(obj,x)
            if (isempty(obj.xhistory))
                obj.xhistory=x';
            else
                obj.xhistory(end+1,:)=x';
            end
        end
        
        function xhistrmsds = getxhistoryrmsdsto(obj,coords)
            for i=1:length(obj.getxhistorycoords()),
                xhistorycoords = obj.getxhistorycoords();
                xhistrmsds(i) = xhistorycoords{i}.registerto(coords);
            end
        end
        
        function xhistcoords = getxhistorycoords(obj)
            if (isempty(obj.xhistcoords)),
                dprintf('converting optimization geometry history (%i iterations), please wait..',size(obj.xhistory,1));
                xhistcoords = cell(size(obj.xhistory,1),1);
                for i=1:size(obj.xhistory,1),
                    xcurr = obj.xhistory(i,:)';
                    xcurrds = data.SimpleDataset([xcurr(1:obj.nfrags) xcurr(2:obj.nfrags+1) obj.geompred(obj.nfrags+2:end)],{'c','cp','z'});            

                    beta = xcurrds.betatransform(obj.sseq,obj.optidata.sclass.classifier);
                    predbond = geom.geometry2bond(xcurr(1:obj.nfrags+1),cos(beta),sin(beta),obj.avglen*ones(obj.nfrags+2,1));

                    xhistcoords{i} = geom.Coords(predbond);
                end
                obj.xhistcoords = xhistcoords;
            else
                xhistcoords = obj.xhistcoords;
            end
        end
        
        function preparebounds(obj,optidata)
            obj.nfrags = length(obj.sgseq);
            % x = [c_1, c_2, ..., c_nfrags, c_nfrags+1, z_1,..,z_nfrags]
            obj.dimx = 2*obj.nfrags + 1; % bond and torsion angles                                                
            obj.ind_c = 1:obj.nfrags+1;
            obj.ind_z = obj.nfrags + 2 : obj.dimx;
            % generate objective    
            covseq = optidata.sgcovmodel(obj.sgseq);
            obj.fevaldata = opti.FEvalData(covseq,obj.aaseq,obj.sseq,obj.sclassifier,obj.options.ccmax);
            if (obj.options.fastobjective)
                obj.fobj = @(x)obj.fevaldata.Vseq(x);                
                %[ nfrags, idx, idx2, means, Ls, individualents ] = opti.seq_pots_fast_prepare( covseq );
                %obj.fobj = @(x) opti.seq_pots_fast(x, nfrags, idx, idx2, means, Ls, individualents);                
            elseif (obj.options.savexhistory)
                obj.fobj = @(x) opti.seq_pots(x,covseq,@obj.savexcallback);                            
            else
                obj.fobj = @(x) opti.seq_pots(x,covseq);
            end
            % generate bounds and starting values
            % for bounds: intersect
            % for X0: average c/cp mus from potential and use t mu
            % TODO: X0 averageing could be done better as we have more information (cov) that we could use

            % LB <= X <= UB
            obj.LB = zeros(obj.dimx,1);
            obj.UB = obj.LB;
            obj.X0 = obj.UB;
            % for c
            for i = obj.ind_c,    

                if (i>1 && i<=obj.nfrags)
                    frag_sgamma_i = obj.sgseq(i-1);
                    frag_sgamma_ip = obj.sgseq(i);
                    % intersect bounds: c_i+1 \in \c(sgamma_i) and c_i \in \c(sgamma_i+1) 
                    obj.LB(i) = max(optidata.sgbounds.cplower(frag_sgamma_i),optidata.sgbounds.clower(frag_sgamma_ip));
                    obj.UB(i) = min(optidata.sgbounds.cpupper(frag_sgamma_i),optidata.sgbounds.cupper(frag_sgamma_ip));

                    % take weighted average of means
                    pot_i = optidata.sgcovmodel(frag_sgamma_i);
                    pot_ip = optidata.sgcovmodel(frag_sgamma_ip);        
                    obj.X0(i) = (pot_i.mean(2)*pot_i.freq + pot_ip.mean(1)*pot_ip.freq) / (pot_ip.freq + pot_i.freq);        
                elseif i==1 % first
                    frag_sgamma_i = obj.sgseq(1);        
                    obj.LB(i) = optidata.sgbounds.clower(frag_sgamma_i);
                    obj.UB(i) = optidata.sgbounds.cupper(frag_sgamma_i);
                    obj.X0(i) = optidata.sgcovmodel(frag_sgamma_i).mean(1);
                else % last
                    frag_sgamma_i = obj.sgseq(end);        
                    obj.LB(i) = optidata.sgbounds.cplower(frag_sgamma_i);
                    obj.UB(i) = optidata.sgbounds.cpupper(frag_sgamma_i);
                    obj.X0(i) = optidata.sgcovmodel(frag_sgamma_i).mean(2);        
                end

                if (obj.LB(i) > obj.UB(i))
                    dprintf('empty bound interval for bond angle c_%i',i)
                end
            end

            % bounds on z
            % TODO: absolute bounds ("both")
            for i = obj.ind_z,    
                frag_sgamma_i = obj.sgseq(i-obj.nfrags-1);
                obj.LB(i) = optidata.sgbounds.zlower(frag_sgamma_i);
                obj.UB(i) = optidata.sgbounds.zupper(frag_sgamma_i);   
                obj.X0(i) = optidata.sgcovmodel(frag_sgamma_i).mean(3);
            end
            
            % overwrite X0 with supplied starting point
            if (obj.options.useX0)
                obj.X0 = obj.options.X0;
            end
            
            % generate linear inequality constraints
            obj.A = [];
            obj.B = [];

            % generate equality constraints
            obj.Aeq = [];
            obj.Beq = [];

            %generate nonlinear (in)equality constraint function
            if (any([obj.options.usehydroconstr, obj.options.usevmaxconstr, obj.options.useccconstr ]))                     
                %[ nfrags, idx, means, Ls, individualents ] = seq_pots_speed_prepare( covseq );
                %obj.NONLINCON = @(x)opti.nonlinconstr_fast(x, seqmaxV, nfrags, idx, idx2, means, Ls, individualents);
                
                obj.fevaldata.seqmaxV = optidata.sgmaxV(obj.sgseq);
                obj.fevaldata.usehydroconstr = obj.options.usehydroconstr;
                obj.fevaldata.usevmaxconstr = obj.options.usevmaxconstr;
                obj.fevaldata.useccconstr = obj.options.useccconstr;
%                 obj.fevaldata.hydropotseq = obj.options.hydro_pot(obj.aaseq);
%                 obj.fevaldata.hydro_chin = obj.options.hydro_chin;
%                 obj.fevaldata.hydro_min = obj.options.hydro_min;
%                 obj.fevaldata.hydro_max = obj.options.hydro_max;
                obj.fevaldata.gradient = obj.options.gradient;
                
                obj.NONLINCON = @(x)obj.fevaldata.nonlinconstr(x);                
            else
                obj.NONLINCON = [];
            end
        end
        
        function obj = run(obj,gettime)
            %[prediction,FVAL,EXITFLAG,OUTPUT,LAMBDA] = opti.runknitro( chainsgclass.cl, optidata );            
            % KTRLINK(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS,knitroOptionsFile)
            if nargin<2,
                gettime = 1;
            end
            obj.havestarted = 1;
            if (gettime), tic; end
            [obj.geompred,obj.FVAL,obj.EXITFLAG,obj.OUTPUT,obj.LAMBDA] =...
                ktrlink(obj.fobj,obj.X0,obj.A,obj.B,obj.Aeq,obj.Beq,obj.LB,obj.UB,obj.NONLINCON,obj.options.ktroptions);
            if (gettime), obj.knitrotime = toc; end
            obj.haverun = 1;
            
            % transform solution
            obj.preddataset = data.SimpleDataset([obj.geompred(1:obj.nfrags) obj.geompred(2:obj.nfrags+1) obj.geompred(obj.nfrags+2:end)],{'c','cp','z'});            
            
            beta = obj.preddataset.betatransform(obj.sseq,obj.sclassifier);
            predbond = geom.geometry2bond(obj.geompred(1:obj.nfrags+1),cos(beta),sin(beta),38*ones(obj.nfrags+2,1));
            
            obj.predcoords = geom.Coords(predbond);
            
        end
        
    end
    
end

