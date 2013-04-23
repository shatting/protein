classdef entropyopt
    %ENTROPYOPT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pot
        data
        cl
        potval_options
    end
    
    methods
        
        function obj = entropyopt(pot, data, cl, potval_options)
           obj.pot = pot;
           obj.data = data;
           obj.cl = cl;
           obj.potval_options = potval_options;                      
           
           startup_lib('cvx');
        end
        
        function [ entout, delta ] = getentropies(obj, method)
            switch(method)
                case 1 
                    [ entout, delta ] = obj.getentropies_method1();
                case 2 
                    [ entout, delta ] = obj.getentropies_method2();
                case 3 
                    [ entout, delta ] = obj.getentropies_method3();
            end
        end
        function [ entout, delta ] = getentropies_method1(obj)
            %GETENTROPIES Summary of this function goes here
            %   Detailed explanation goes here
           Vgd = suffclass.entropyopt.vgammadelta( obj.pot, obj.data, obj.cl, obj.potval_options );
           
            
            ncl = size(Vgd,1);

            % empty classes. cvx doesnt like inf in constraints
            idxninf = max(triu(Vgd,1),[],2) ~= inf;
            nclninf = sum(idxninf);

            Vgd = Vgd(idxninf,idxninf);

            cvx_begin
                variables delta ent(nclninf);
                maximize(delta);
                subject to
                    for g=1:nclninf,
                        for d=1:nclninf,
                         if d~=g,
                            ent(g) - ent(d) + delta <= Vgd(g,d);  
                         end;
                        end
                    end
                    ent(nclninf) == 0;                    
            cvx_end

            % only one loop try
            % triag = triu(ones(nclninf),1);
            % lin = zeros(nclninf);
            % lin(triag==1) = 1:(nclninf)*(nclninf-1)/2;
            % 
            % S = sparse(lin);
            % [g,d,s] = find(S);
            % 
            % nineq = length(g);
            % 
            % cvx_begin
            %     variables delta ent(nclninf);
            %     maximize delta;
            %     subject to
            %         for i=1:nineq,
            %             ent(g(i)) - ent(d(i)) + delta <= Vgd(g(i),d(i));
            %         end
            %         ent <= 10;
            %         ent(1) == 0;
            % cvx_end
            %----

            entout = zeros(ncl,1) + inf;
            entout(idxninf) = ent;
            
            if (delta < 0)
                dprintf('method 1 found no optimal solution: delta = %g', delta);
                entout
                dprintf('trying method 2.');
                [ entout, delta ] = obj.getentropies_method2;
            end

        end
        
        function [ entout, delta ] = getentropies_method2(obj)
            
            Vgd = suffclass.entropyopt.vgammadelta( obj.pot, obj.data, obj.cl, obj.potval_options );
            
            ncl = size(Vgd,1);

            % empty classes. cvx doesnt like inf in constraints
            idxninf = max(triu(Vgd,1),[],2) ~= inf;
            nclninf = sum(idxninf);

            Vgd = Vgd(idxninf,idxninf);
            
            epsilon = zeros(1,nclninf);
            for g = 1:nclninf,
                epsilon(g) = 1e-8*median(Vgd(g,:));
            end
            
            cvx_begin
                variables delta(nclninf) ent(nclninf);
                maximize(sum(delta));
                subject to
                    for g=1:nclninf,
                        for d=1:nclninf,
                             if d~=g,
                                ent(g) - ent(d) + delta(g) <= Vgd(g,d);  
                             end;
                        end
                        delta(g) <= epsilon(g);
                    end
                    ent(1) == 0;                         
            cvx_end

            
            entout = zeros(ncl,1) + inf;
            entout(idxninf) = ent;
        end
        
        function [ entout, delta ] = getentropies_method3(obj)
            
            datai = obj.data;
            cli = obj.cl;
            delta = -inf;            
            pmisclassified = 1;
            ncl = obj.pot.suff.ncl;
            it = 1;
            deltas = [];
                        
            while delta < 0 && pmisclassified > 0.01
                deltaold = delta;
                Vgd = suffclass.entropyopt.vgammadelta( obj.pot, datai , cli, obj.potval_options );

                % empty classes. cvx doesnt like inf in constraints
                idxninf = max(triu(Vgd,1),[],2) ~= inf;
                nclninf = sum(idxninf);

                Vgd = Vgd(idxninf,idxninf);

                cvx_begin
                    cvx_precision low
                    cvx_quiet(true)
                    variables delta ent(nclninf);                                
                    maximize(delta);
                    subject to
                        for g=1:nclninf,
                            for d=1:nclninf,
                             if d~=g,
                                ent(g) - ent(d) + delta <= Vgd(g,d);  
                             end;
                            end
                        end
                        ent(1) == 0;                    
                cvx_end

                entout = zeros(ncl,1) + inf;
                entout(idxninf) = ent;
                
                potvals = suffclass.potential_values(obj.pot,obj.potval_options);
                potvals.calcvals(datai);
                V = potvals.V;
                Veffective = zeros(size(V,1),1);
                entvals = zeros(size(V));
                for i=1:ncl,
                    Veffective(cli==i) = V(cli==i,i) + ent(i);
                    entvals(:,i) = ent(i);
                end                
                deltal = min(V + entvals,[],2) - Veffective;
                %discard worst                
                idxmisclassified = find(deltal<0);                
                nmisclassified = length(idxmisclassified);
                [a,b] = sort(deltal(deltal<0));
                
                cutidx = round(0.05*nmisclassified);
                if (cutidx == 1)
                    break;
                end                
                idxtokeep = cutidx:nmisclassified;
                
                
                %now transpose idx back into datai, cli
                oidxtokeep = b(idxtokeep);
                                
                idxkeep = deltal==0; % keep all correctly classifieds
                idxkeep(idxmisclassified(oidxtokeep)) = 1; % ..and the ones we decided to keep
                
                % apply filtering
                %dprintf('Iteration %i done.\nHad %i(%.2f%%) misclassified.',it,nmisclassified,100*nmisclassified/length(deltal));
                %dprintf('keeping %i of %i',sum(idxkeep),length(idxkeep));
                datai = datai(idxkeep,:);
                cli = cli(idxkeep);
                deltas = [deltas delta];                        
                it = it +1;
                pmisclassified = nmisclassified/length(deltal);
            end
        end
        
    end
    
    methods (Static)
        
        Vgd = vgammadelta( pot, data, cl, options )
            
    end
    
end

