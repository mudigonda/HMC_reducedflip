function dEdX = dEdX_gauss( X, J, mu )
% return the gradient with respect to X for the energy of each sample
% (column vector) in X for a Gaussian with inverse covariance matrix
% (coupling matrix) J

    dEdX = 0;

    px = 0;

    for ii = 1:length(J)
        Jl = J{ii};
        Xl = bsxfun( @plus, X, repmat(-mu{ii},size(X,2),1)');
        px = px + exp(-0.5*sum( Xl.*(Jl*Xl), 1 ));
    end

    E = -log(px);
    
    for ii = 1:length(J)
        Jl = J{ii};
        Xl = bsxfun( @plus, X, repmat(-mu{ii},size(X,2),1)');
        dEdX = dEdX + exp(-0.5*sum( Xl.*(Jl*Xl), 1 ))*(0.5*Jl*Xl + 0.5*Jl'*Xl);
    end
    
    dEdX = dEdX/px;
    
    if ~isfinite(E)
        dEdX(:) = 0;
    end
    
%    if sum(~isfinite(dEdX) ) > 0
%        keyboard
%    end