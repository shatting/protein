classdef FEvalData < handle
    %FEVALDATA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
         nfrags
         idx
         idx2
         means
         Ls
         e_i
         
         % constraint options
         usevmaxconstr
         seqmaxV
         
         usehydroconstr
         useccconstr        
         sseq
         aaseq
         sclassifier
         ccmax
         cosbetabar
         sinbetabar
         hydropotseq
         hydro_chin
         hydro_min
         hydro_max
         
         gradient
                  
         % lastx stuff
         lastxenable = 0;
         
         lastx
         lastxfval                  
         lastxpots

         lastxgradc1 % L_{:1}^(j)^T L^(j) (x_j-mu_j)
         lastxgradc2 % L_{:2}^(j)^T L^(j) (x_j-mu_j)
         lastxgradz  % L_{:3}^(j)^T L^(j) (x_j-mu_j)
         lastxgrad
    end
    
    methods
        function obj = FEvalData(potseq,aaseq,sseq,sclassifier,ccmax)
            obj.aaseq = aaseq;
            obj.sseq = sseq;
            obj.sclassifier = sclassifier;            
            obj.cosbetabar = cos(sclassifier.gettorsioncenters)';
            obj.sinbetabar = sin(sclassifier.gettorsioncenters)';
            obj.ccmax = ccmax;
            [  obj.nfrags, obj.idx, obj.idx2, obj.means, obj.Ls, obj.e_i ] = obj.initialize( potseq );
        end
                
        % function left here for custom evaluation
        function [  nfrags, idx, idx2, means, Ls, e_i ] = initialize(obj, potseq )            
            
            nfrags = length(potseq);
            % want: [c_1,c_2,z_1 ; c_2,c_3,z_2 ; ..]'
            idx = [(1:nfrags)',(2:nfrags+1)',(nfrags+2:2*nfrags+1)']';
            e_i = [potseq.ent];
            means = [potseq.mean];
            if 0, % fast option
                Ls = zeros(3,3,nfrags);
                for i=1:nfrags,
                    Ls(:,:,i) = potseq(i).L;
                end
                idx2 = [];
            else % turbo option
                idx2 = 1:nfrags;
                idx2 = idx2([idx2;idx2;idx2]);    
                idx2 = reshape(idx2,3*nfrags,1);

                Ls = zeros(3*nfrags,3);
                for i=1:nfrags,
                    Ls(3*(i-1)+1:3*i,:) = potseq(i).L;
                end            
            end
        end
                
        function [ fval, grad, V_i, gradc1, gradc2, gradz ] = Vseq(obj,x,givegrad)
            % same as seq_pots/_fast
            %  (0.938 secs->0.639 secs for all chains in db (fast))
            %  (0.938 secs->0.134 secs for all chains in db (turbo))         
            %
            % INPUT:
            %       x           [c_1,..,c_(n-2),z_1,..,z_(n-3)]'
            %              
            % OUTPUT:
            %       fval        sum_i(V_{gamma,s}(x_{i:i+2} + ent_{gamma,s})
            %       V_i 

%             if (obj.lastxenable && ~isempty(obj.lastx) && all(x==obj.lastx))
%                 fval = obj.lastxfval;
%                 grad = obj.lastxgrad;                
%                 if (nargout>2)       
%                     V_i = obj.lastxpots;
%                     gradc1 = obj.lastxgradc1;
%                     gradc2= obj.lastxgradc2;
%                     gradz= obj.lastxgradz;
%                 end
%                 return
%             end
            
            xv = x(obj.idx);
            xminusmu = (xv - obj.means)';
            Lxminusmu = obj.Ls.*xminusmu(obj.idx2,:);            
            Lxminusmu = sum(Lxminusmu,2); % = [L^(1) (x_1-mu_1);L^(2) (x_2-mu_2);...]            
            V_i = reshape(Lxminusmu,3,obj.nfrags);
            V_i = sum(V_i.^2,1);

            fval = sum(V_i + obj.e_i);                               
            
            if (obj.lastxenable || nargout > 1 && (nargin < 3 || givegrad))
                %compute gradient
                gradc1 = [obj.Ls(:,1).*Lxminusmu;0;0;0]; % nfrags-2 x 1
                gradc2 = [0;0;0;obj.Ls(:,2).*Lxminusmu]; % nfrags-2 x 1
                gradz = obj.Ls(:,3).*Lxminusmu; % nfrags-3 x 1
                grad = reshape([gradc1 + gradc2;gradz],3,2*obj.nfrags+1);
                grad = 2*sum(grad,1)';
            else
                grad = [];
            end
            
%             if (obj.lastxenable)
%                 obj.lastx = x;
%                 obj.lastxfval = fval;
%                 obj.lastxpots = V_i;
%                 if (~isempty(grad))
%                     obj.lastxgrad = grad;
%                     obj.lastxgradc1 = gradc1;
%                     obj.lastxgradc2 = gradc2;
%                     obj.lastxgradz = gradz;
%                 end
%             end
        end

        function [ c, ceq, cgrad, ceqgrad ] = nonlinconstr(obj, x)            
                        
            if (nargout > 2) % gradients requested
                [ ~, ~, V_i, gradc1, gradc2, gradz ] = obj.Vseq(x,1);                                    
                
                c = (V_i + obj.e_i) - obj.seqmaxV;
                
                gradc1 = reshape(gradc1(1:end-3),3,obj.nfrags);
                gradc2 = reshape(gradc2(4:end),3,obj.nfrags);
                gradz = reshape(gradz,3,obj.nfrags);
                
                gradc1 = sum(gradc1,1);
                gradc2 = sum(gradc2,1);
                gradz = sum(gradz,1);                
                
                nf=obj.nfrags;
                cgrad = zeros(2*nf+1,nf)';
                idx3 = 1:nf+1:nf^2;
                cgrad(idx3) = gradc1;
                cgrad(idx3+nf) = gradc2;
                cgrad(idx3+nf*(nf+1)) = gradz;
                cgrad = 2*cgrad';
                
                % get coords
                if (obj.useccconstr || obj.usehydroconstr)
                    xgrad = gradientinit(x);
                    [t,tp] = geom.ttptransform_rational(xgrad(obj.nfrags+2:end),obj.sseq,obj.cosbetabar,obj.sinbetabar);                    
                    bond = geom.geometry2bond(xgrad(1:obj.nfrags+1),t,tp,38*ones(obj.nfrags+2,1),strcmp(obj.gradient,'intlab'));
                    coords = geom.bond2coords(bond);
                end
                
                if (obj.useccconstr)
                    cc = opti.constr_cc(coords);
                    c = [c,cc-obj.ccmax];
                    cgrad = [cgrad, full(c(1,end-2).dx')];
                end

                if (obj.usehydroconstr)
                    vhydro = opti.constr_hydro(coords,obj.hydropotseq,obj.hydro_chin);
                    c = [c,vhydro - obj.hydro_max,-vhydro + obj.hydro_min];                
                    cgrad = [cgrad, full(c(1,end-1).dx'),full(c(1,end).dx')];
                end
                                
                c=c.x;
                
                ceqgrad = [];
            else
                [~, ~, seqVraw] = obj.Vseq(x,0); % we want pots not grad
                c = (seqVraw + obj.e_i) - obj.seqmaxV; 
                
                % get coords
                if (obj.useccconstr || obj.usehydroconstr)
                    ci = x(1:obj.nfrags+1);
                    %preddataset = data.SimpleDataset([ci(1:end-1) ci(2:end) x(obj.nfrags+2:end)],{'c','cp','z'});                        
                    %[t,tp] = preddataset.ttptransform_rational(obj.sseq,obj.sclassifier);
                    [t,tp] = geom.ttptransform_rational(x(obj.nfrags+2:end),obj.sseq,obj.cosbetabar,obj.sinbetabar);
                    bond = geom.geometry2bond(ci,t,tp,38*ones(obj.nfrags+2,1));
                    coords = geom.bond2coords(bond);                                                
                end
                
                if obj.useccconstr
                    cc = opti.constr_cc(coords);
                    c = [c,cc-obj.ccmax];
                end

                if (obj.usehydroconstr)
                    vhydro = opti.constr_hydro(coords,obj.hydropotseq,obj.hydro_chin);
                    c = [c,vhydro - obj.hydro_max,-vhydro + obj.hydro_min];                
                end

            end           
            
            if (0 && any(c)>0)
                dprintf('quadratic constraint violated!');
            end
            ceq = [];            

        end
        
    end
    
end

