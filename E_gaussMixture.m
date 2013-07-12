function E = E_gaussMixture( X, J, mu )
% return the energy for each sample (column vector) in X for a Gaussian with
% inverse covariance matrix (coupling matrix) J
   
    px = 0;

    for ii = 1:length(J)
        Jl = J{ii};
        Xl = bsxfun( @plus, X, -mu{ii} );
        px = px + exp(-0.5*sum( Xl.*(Jl*Xl), 1 ));
    end

    E = -log(px);