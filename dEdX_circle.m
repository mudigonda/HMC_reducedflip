function dEdX = dEdX_circle( X, lambda )
% return the gradient with respect to X for the energy of each sample
% (column vector) in X for a Gaussian with inverse covariance matrix
% (coupling matrix) J

	scl = 0.1;
	X = X*scl;

    u = sum(X.^2, 1 );
    r = sqrt( u );
    %E = lambda * log(r).^2;

    dEdr = 2*lambda*log(r)./r;
    dEdu = (1/2)./r.*dEdr;
    dEdX = bsxfun( @times, dEdu, 2*X );

    dEdX = dEdX*scl;
