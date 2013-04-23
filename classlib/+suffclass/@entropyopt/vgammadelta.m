function [ Vgd ] = vgammadelta( pot, data, cl, options )
%VGAMMADELTA get V_{\gamma\delta} := min_l{V_\delta(y_l) - v_\gamma(y_l) | \gamma = \gamma_l}

% get raw potential values, without entropy or priors, in V
potvals = suffclass.potential_values(pot,options);
potvals.calcvals(data);
V = potvals.V;

Vgd = zeros(pot.suff.ncl) + inf;
Vgd(logical(eye(pot.suff.ncl))) = 0;

for gamma=1:pot.suff.ncl,
	indgamma = cl==gamma;
    
    if (sum(indgamma)==0), continue; end
    
    for delta = setdiff(1:pot.suff.ncl,gamma),
                
        Vgd(gamma,delta) = min(V(indgamma,delta) - V(indgamma,gamma));
        
    end
    
end


end
