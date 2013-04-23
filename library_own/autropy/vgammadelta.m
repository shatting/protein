function [ Vgd ] = vgammadelta( pot, data, cl )
%VGAMMADELTA get V_{\gamma\delta} := min_l{V_\delta(y_l) - v_\gamma(y_l) | \gamma = \gamma_l}

% get raw potential values, without entropy or priors, in V
[ cpred, potmin, V ] = suff_potvals( pot, data );

Vgd = zeros(pot.ncl) + inf;

for gamma=1:pot.ncl,
	indcl = cl==gamma;
    
    if (sum(indcl)==0), continue; end
    
    for delta = 1:pot.ncl,
        
        diff = V(indcl,delta) - V(indcl,gamma);
        Vgd(gamma,delta) = min(diff);
        
    end
    
end


end
