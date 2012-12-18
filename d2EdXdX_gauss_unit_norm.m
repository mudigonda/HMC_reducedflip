function d2E = d2EdXdX_gauss_unit_norm( X, varargin )
    
    szb = size(X,2);
    szd = size(X,1);
    
    d2E = bsxfun( @times, reshape(eye(szd), [szd,szd,1] ), reshape( ones( szb,1), [1,1,szb] ) );
    