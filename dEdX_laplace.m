function dE = dEdX_laplace( X, theta )
% product of laplace experts
ndims = size(X, 1);
nbatch = size(X, 2);
nexperts = prod(size(theta)) / (ndims);
W = reshape( theta, [nexperts, ndims] );
ff = W * X;
dE = W' * sign(ff);
end