% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)
% copyright 2010 Jascha Sohl-Dickstein

function E = E_gauss_unit_norm( x, varargin )

%E = 0.5 * x' * x;
    E = 0.5 * sum( x.^2, 1 );
    