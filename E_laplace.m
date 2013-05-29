function E = E_laplace( X, theta )
% product of laplace experts
ndims = size(X, 1);
nbatch = size(X, 2);
nexperts = prod(size(theta)) / (ndims);
W = reshape( theta, [nexperts, ndims] );
ff = W * X;
E = sum( abs(ff), 1 );
end