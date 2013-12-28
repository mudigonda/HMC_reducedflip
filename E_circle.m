function E = E_circle( X, lambda )
% return the energy for each sample (column vector) in X for a Gaussian with
% inverse covariance matrix (coupling matrix) J

	scl = 0.1;
	X = X*scl;

    r = sqrt( sum(X.^2, 1 ) );
    E = lambda * log(r).^2;
